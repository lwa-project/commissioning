#!/usr/bin/env python3
    
import os
import sys
import aipy
import copy
import numpy as np
import argparse
from calendar import timegm
from datetime import datetime

from lsl import astro
from lsl.common import stations
from lsl.statistics.robust import *
from lsl.correlator import uvutils
from lsl.writer.fitsidi import NUMERIC_STOKES
from lsl.misc import parser as aph

from lsl.skymap import SkyMapLFSM, ProjectedSkyMap
from lsl.imaging import utils, selfcal, overlay
from lsl.sim import vis as simVis

import logging
from lsl.logger import enable_file_logging, set_log_level
set_log_level(logging.DEBUG)

import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter


def main(args):
    if args.log:
        enable_file_logging(args.log)
        
    filename = args.filename
    
    idi = utils.CorrelatedData(filename)
    aa = idi.get_antennaarray()
    lo = idi.get_observer()
    lo.date = idi.date_obs.strftime("%Y/%m/%d %H:%M:%S")
    jd = lo.date + astro.DJD_OFFSET
    lst = str(lo.sidereal_time())
    
    nStand = len(idi.stands)
    nchan = len(idi.freq)
    freq = idi.freq
    
    print("Raw Stand Count: %i" % nStand)
    print("Final Baseline Count: %i" % (nStand*(nStand-1)//2,))
    print("Spectra Coverage: %.3f to %.3f MHz in %i channels (%.2f kHz/channel)" % (freq[0]/1e6, freq[-1]/1e6, nchan, (freq[-1] - freq[0])/1e3/nchan))
    print("Polarization Products: %i starting with %i" % (len(idi.pols), idi.pols[0]))
    print("JD: %.3f" % jd)
    
    # Pull out something reasonable
    toWork = np.where((freq>=args.lower) & (freq<=args.upper))[0]
    
    print("Reading in FITS IDI data")
    nSets = idi.total_baseline_count // (nStand*(nStand+1)//2)
    for set in range(1, nSets+1):
        print("Set #%i of %i" % (set, nSets))
        fullDict = idi.get_data_set(set, include_auto=True)
        autoDict = fullDict.get_uv_range(max_uv=0.001)
        
        if args.filter_autos:
            # Filter out dead or high power antennas based on the autos
            auto_pwr = {}
            for i,(ax,ay) in enumerate(zip(autoDict.XX.data, autoDict.YY.data)):
                auto_pwr[i] = [np.abs(ax[toWork]), np.abs(ay[toWork])]
            pwr = np.array(list(auto_pwr.values()))
            pwr = np.log10(np.clip(pwr, 1e-10, np.inf))
            med = np.median(pwr, axis=0)
            mad = np.median(np.abs(pwr-med), axis=0)
            flg = (np.abs(pwr[:,0,:]-med[0,:]) > 5*mad[0,:]*1.4826).mean(axis=1) * 0.5 \
                + (np.abs(pwr[:,1,:]-med[1,:]) > 5*mad[1,:]*1.4826).mean(axis=1) * 0.5
            bad = np.where(flg > 0.3)[0]
            if len(bad):
                print(f"Flagging {len(bad)} antennas based off auto-correlations:")
                for b in bad:
                    print(f"  Stand {idi.stands[autoDict.baselines[b][0]]} with {flg[b]:.0%} of channels-pols flagged")
                bad_ants = [autoDict.baselines[b][0] for b in bad]
                fullDict = fullDict.get_antenna_subset(exclude=bad_ants, indicies=True)
                
        # Downselect
        dataDict = fullDict.get_uv_range(min_uv=0.001, max_uv=args.max_uv_dist)
        dataDict.sort()
        
        # Gather up the polarizations and baselines
        pols = dataDict.pols
        bls = dataDict.baselines
        print("The reduced list has %i baselines and %i channels" % (len(bls), len(toWork)))
        
        # Build a list of unique JDs for the data
        jdList = [dataDict.jd]
        
        # Build the simulated visibilities
        print("Building Model")
        sm = SkyMapLFSM(freq_MHz=freq.mean()/1e6)
        pm = ProjectedSkyMap(sm, aa.lat*180/np.pi, aa.lon*180/np.pi, jdList[0])
        pm.compute_visible_power()
        SOURCES = {}
        for i in range(pm.visibleRa.size):
            SOURCES[f"pnt{i+1}"] = aipy.amp.RadioFixedBody(pm.visibleRa[i]*np.pi/180,
                                                           pm.visibleDec[i]*np.pi/180,
                                                           jys=pm.visibleNormalizedPower[i],
                                                           index=-2.3)
        print(f"Using a source catalog with {len(SOURCES)} entries")
        set_log_level(logging.INFO)     # Too much info with DEBUG in the next step
        simDict = simVis.build_sim_data(aa, SOURCES, jd=[jdList[0],], pols=pols, baselines=bls)
        set_log_level(logging.DEBUG)
        
        print("Running self cal.")
        simDict.sort()
        dataDict.sort()
        fixedDataXX, delaysXX, convXX = selfcal.delay_only(aa, dataDict, simDict, toWork, 'XX',
                                                           ref_ant=args.reference,
                                                           max_iter=args.max_iterations,
                                                           delay_cutoff=args.delay_cutoff,
                                                           return_convergence=True)
        fixedDataYY, delaysYY, convYY = selfcal.delay_only(aa, dataDict, simDict, toWork, 'YY',
                                                           ref_ant=args.reference,
                                                           max_iter=args.max_iterations,
                                                           delay_cutoff=args.delay_cutoff,
                                                           return_convergence=True)
        
        print("    Saving results")
        outname = os.path.split(filename)[1]
        outname = os.path.splitext(outname)[0]
        outname = "%s.sc" % outname
        with open(outname, 'w') as fh:
            fh.write("################################\n")
            fh.write("#                              #\n")
            fh.write("# Settings:                    #\n")
            fh.write("#  Method: diffuse             #\n")
            fh.write(f"#  Ref Ant: {args.reference:4d}               #\n")
            fh.write(f"#  Lower: {args.lower/1e6:5.1f} MHz            #\n")
            fh.write(f"#  Upper: {args.upper/1e6:5.1f} MHz            #\n")
            fh.write(f"#  Filtering: {str(args.filter_autos):5s}            #\n")
            fh.write(f"#  Max (u,v): {args.max_uv_dist:4.1f} lambda      #\n")
            fh.write(f"#  Max Iters: {args.max_iterations:3d}              #\n")
            fh.write(f"#  Delay Cutoff: {args.delay_cutoff:4.2f} ns       #\n")
            fh.write(f"#  Inv Eps: {args.inv_epsilon:5.2f}              #\n")
            fh.write("#                              #\n")
            fh.write(f"#  Converged XX: {str(convXX):5s}         #\n")
            fh.write(f"#  Converged YY: {str(convYY):5s}         #\n")
            fh.write("#                              #\n")
            fh.write("################################\n")
            fh.write("#                              #\n")
            fh.write("# Columns:                     #\n")
            fh.write("# 1) Stand number              #\n")
            fh.write("# 2) X pol. amplitude          #\n")
            fh.write("# 3) X pol. delay (ns)         #\n")
            fh.write("# 4) Y pol. amplitude          #\n")
            fh.write("# 5) Y pol. delay (ns)         #\n")
            fh.write("#                              #\n")
            fh.write("################################\n")
            if args.filter_autos:
                fh.write("#                              #\n")
                fh.write("# Filtered Stands:             #\n")
                for b in bad:
                    fh.write(f"#  Stand {idi.stands[autoDict.baselines[b][0]]:3d}                   #\n")
                fh.write("#                              #\n")
                fh.write("################################\n")
            for i in range(delaysXX.size):
                fh.write("%3i  %.6g  %.6g  %.6g  %.6g\n" % (idi.stands[i], 1.0, delaysXX[i], 1.0, delaysYY[i]))
                
        # Build up the images for each polarization
        if args.plot:
            fixedFullXX = simVis.scale_data(fullDict, delaysXX*0+1, delaysXX)
            fixedFullYY = simVis.scale_data(fullDict, delaysYY*0+1, delaysYY)
            
            print("    Gridding")
            toWork = np.where((freq>=30e6) & (freq<=82e6))[0]
            try:
                imgXX = utils.build_gridded_image(fullDict, size=80, res=0.5, pol='XX', chan=toWork)
            except:
                imgXX = None
                
            try:
                imgYY = utils.build_gridded_image(fullDict, size=80, res=0.5, pol='YY', chan=toWork)
            except:
                imgYY = None
                
            try:
                simgXX = utils.build_gridded_image(simDict, size=80, res=0.5, pol='XX', chan=toWork)
            except:
                simgXX = None
            try:
                simgYY = utils.build_gridded_image(simDict, size=80, res=0.5, pol='YY', chan=toWork)
            except:
                simgYY = None
                
            try:
                fimgXX = utils.build_gridded_image(fixedFullXX, size=80, res=0.5, pol='XX', chan=toWork)
            except:
                fimgXX = None
            try:
                fimgYY = utils.build_gridded_image(fixedFullYY, size=80, res=0.5, pol='YY', chan=toWork)
            except:
                fimgYY = None
                
            # Plots
            print("    Plotting")
            fig = plt.figure()
            ax1 = fig.add_subplot(3, 2, 1)
            ax2 = fig.add_subplot(3, 2, 2)
            ax3 = fig.add_subplot(3, 2, 3)
            ax4 = fig.add_subplot(3, 2, 4)
            ax5 = fig.add_subplot(3, 2, 5)
            ax6 = fig.add_subplot(3, 2, 6)
            for ax, img, pol in zip([ax1, ax2, ax3, ax4, ax5, ax6], [imgXX, imgYY, simgXX, simgYY, fimgXX, fimgYY], ['XX', 'YY', 'simXX', 'simYY', 'scalXX', 'scalYY']):
                # Skip missing images
                if img is None:	
                    ax.text(0.5, 0.5, 'Not found in file', color='black', size=12, horizontalalignment='center')
                    
                    ax.xaxis.set_major_formatter( NullFormatter() )
                    ax.yaxis.set_major_formatter( NullFormatter() )
                    
                    ax.set_title("%s @ %s LST" % (pol, lst))
                    continue
                    
                # Display the image and label with the polarization/LST
                out = img.image(center=(80,80))
                print(pol, out.min(), out.max())
                #if pol == 'scalXX':
                    #out = np.rot90(out)
                    #out = np.rot90(out)
                cb = ax.imshow(out, extent=(1,-1,-1,1), origin='lower', 
                        vmin=img.image().min(), vmax=img.image().max())
                fig.colorbar(cb, ax=ax)
                ax.set_title("%s @ %s LST" % (pol, lst))
                
                # Turn off tick marks
                ax.xaxis.set_major_formatter( NullFormatter() )
                ax.yaxis.set_major_formatter( NullFormatter() )
                
                # Compute the positions of major sources and label the images
                compSrc = {}
                ax.plot(0, 0, marker='+', markersize=10, markeredgecolor='w')
                for name,src in simVis.SOURCES.items():
                    src.compute(aa)
                    top = src.get_crds(crdsys='top', ncrd=3)
                    az, alt = aipy.coord.top2azalt(top)
                    compSrc[name] = [az, alt]
                    if alt <= 0:
                        continue
                    ax.plot(top[0], top[1], marker='x', markerfacecolor='None', markeredgecolor='w', 
                            linewidth=10.0, markersize=10)
                    ax.text(top[0], top[1], name, color='white', size=12)
                    
                # Add lines of constant RA and dec.
                overlay.horizon(ax, aa)
                overlay.graticule_radec(ax, aa)
                
            plt.show()
            
    print("...Done")


if __name__ == "__main__":
    np.seterr(all='ignore')
    
    parser = argparse.ArgumentParser(
        description="self-calibrate a TBX FITS IDI file using diffuse emission",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
    parser.add_argument('filename', type=str, 
                        help='filename to calibrate')
    parser.add_argument('-r', '--reference', type=aph.positive_int, default=1,
                        help='reference stand to use')
    parser.add_argument('-l', '--lower', type=aph.positive_float, default=35.0,
                        help='lowest frequency to consider in MHz')
    parser.add_argument('-u', '--upper', type=aph.positive_float, default=85.0,
                        help='highest frequency to consider in MHz')
    parser.add_argument('-f', '--filter-autos', action='store_true',
                        help='filter based on the auto-correlation to remove data and weird stands')
    parser.add_argument('-m', '--max-uv-dist', type=aph.positive_float, default=4.0,
                        help='maximimum baseline (u,v) length to use in wavelengths')
    parser.add_argument('-i', '--max-iterations', type=aph.positive_int, default=60,
                        help="maximum number of self-cal iterations")
    parser.add_argument('-d', '--delay-cutoff', type=aph.positive_float, default=0.2,
                        help="delay cutoff in ns for the self-cal convergence threshold")
    parser.add_argument('-e', '--inv-epsilon', type=aph.positive_or_zero_float, default=0.0,
                        help="Tikhonov regularization strength, zero disables")
    parser.add_argument('-p', '--plot', action='store_true',
                        help='plot the results at the end')
    parser.add_argument('--log', type=str,
                        help='save LSL logging output to the specified filename')
    args = parser.parse_args()
    args.lower *= 1e6
    args.upper *= 1e6
    main(args)
