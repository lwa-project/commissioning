#!/usr/bin/env python3

import os
import sys
import aipy
import copy
import pytz
import numpy
import argparse
from calendar import timegm
from datetime import datetime

from lsl import astro
from lsl.common import stations
from lsl.statistics.robust import *
from lsl.correlator import uvutils
from lsl.writer.fitsidi import NUMERIC_STOKES
from lsl.misc import parser as aph

from lsl.imaging import utils, selfcal
from lsl.sim import vis as simVis

import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter


MST = pytz.timezone('US/Mountain')
UTC = pytz.UTC


def graticle(ax, lst, lat, label=True):
    """
    For a matplotlib axis instance showing an image of the sky, plot lines of
    constant declination and RA.  Declinations are spaced at 20 degree intervals
    and RAs are spaced at 2 hour intervals.
    
    .. note::
        LST and latitude values should be passed as radians.  This is the default
        for lwa1.get_observer.sidereal_time() and lwa1.get_observer().lat.
    """
    
    # Lines of constant declination first
    decs = range(-80, 90, 20)
    ras = numpy.linspace(0, 360, 800)
    
    x = numpy.zeros(ras.size)
    x = numpy.ma.array(x, mask=numpy.zeros(ras.size))
    y = numpy.zeros(ras.size)
    y = numpy.ma.array(y, mask=numpy.zeros(ras.size))
    
    for dec in decs:
        x *= 0
        y *= 0
        
        # Loop over RA to compute the topocentric coordinates (used by the image) for
        # the lines.  Also, figure out the elevation for each point on the line so
        # we can mask those below the horizon
        for i,ra in enumerate(ras):
            eq = aipy.coord.radec2eq((-lst + ra*numpy.pi/180,dec*numpy.pi/180))
            xyz = numpy.dot(aipy.coord.eq2top_m(0, lat), eq)
            az,alt = aipy.coord.top2azalt(xyz)
            
            x[i] = xyz[0]
            y[i] = xyz[1]
            if alt <= 0:
                x.mask[i] = 1
                y.mask[i] = 1
            else:
                x.mask[i] = 0
                y.mask[i] = 0
                
        ax.plot(x, y, color='white', alpha=0.75)
        
        eq = aipy.coord.radec2eq((-lst + lst,(dec+5)*numpy.pi/180))
        xyz = numpy.dot(aipy.coord.eq2top_m(0, lat), eq)
        az,alt = aipy.coord.top2azalt(xyz)
        
        if alt > 15*numpy.pi/180 and label:
            ax.text(xyz[0], xyz[1], '%+i$^\circ$' % dec, color='white')
            
    # Lines of constant RA
    decs = numpy.linspace(-80, 80, 400)
    ras = range(0,360,30)
    
    x = numpy.zeros(decs.size)
    x = numpy.ma.array(x, mask=numpy.zeros(decs.size))
    y = numpy.zeros(decs.size)
    y = numpy.ma.array(y, mask=numpy.zeros(decs.size))
    
    for ra in ras:
        x *= 0
        y *= 0
        
        # Loop over dec. to compute the topocentric coordinates (used by the image) for
        # the lines.  Also, figure out the elevation for each point on the line so
        # we can mask those below the horizon
        for i,dec in enumerate(decs):
            eq = aipy.coord.radec2eq((-lst + ra*numpy.pi/180,dec*numpy.pi/180))
            xyz = numpy.dot(aipy.coord.eq2top_m(0, lat), eq)
            az,alt = aipy.coord.top2azalt(xyz)
            
            x[i] = xyz[0]
            y[i] = xyz[1]
            if alt <= 0:
                x.mask[i] = 1
                y.mask[i] = 1
            else:
                x.mask[i] = 0
                y.mask[i] = 0
                
        ax.plot(x, y, color='white', alpha=0.75)
        
        eq = aipy.coord.radec2eq((-lst + ra*numpy.pi/180,0))
        xyz = numpy.dot(aipy.coord.eq2top_m(0, lat), eq)
        az,alt = aipy.coord.top2azalt(xyz)
        
        if alt > 20*numpy.pi/180 and label:
            ax.text(xyz[0], xyz[1], '%i$^h$' % (ra/15,), color='white')


def main(args):
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
    toWork = numpy.where((freq>=args.lower) & (freq<=args.upper))[0]
    
    print("Reading in FITS IDI data")
    nSets = idi.total_baseline_count // (nStand*(nStand+1)//2)
    for set in range(1, nSets+1):
        print("Set #%i of %i" % (set, nSets))
        fullDict = idi.get_data_set(set)
        if args.min_uv_dist > 0.0:
            dataDict = fullDict.get_uv_range(min_uv=args.min_uv_dist)
        else:
            dataDict = fullDict
        dataDict.sort()
        
        # Gather up the polarizations and baselines
        pols = dataDict.pols
        bls = dataDict.baselines
        print("The reduced list has %i baselines and %i channels" % (len(bls), len(toWork)))
        
        # Build a list of unique JDs for the data
        jdList = [dataDict.jd]
         
        # Build the simulated visibilities
        print("Building Model")
        simVis.SOURCES['Sun']._jys *= args.sun_factor
        simDict = simVis.build_sim_data(aa, simVis.SOURCES, jd=[jdList[0],], pols=pols, baselines=bls)
        
        print("Running self cal.")
        simDict.sort()
        dataDict.sort()
        fixedDataXX, delaysXX = selfcal.delay_only(aa, dataDict, simDict, toWork, 'XX',
                                                   ref_ant=args.reference,
                                                   max_iter=args.max_iterations,
                                                   delay_cutoff=args.delay_cutoff)
        fixedDataYY, delaysYY = selfcal.delay_only(aa, dataDict, simDict, toWork, 'YY',
                                                   ref_ant=args.reference,
                                                   max_iter=args.max_iterations,
                                                   delay_cutoff=args.delay_cutoff)
        fixedFullXX = simVis.scale_data(fullDict, delaysXX*0+1, delaysXX)
        fixedFullYY = simVis.scale_data(fullDict, delaysYY*0+1, delaysYY)
        
        print("    Saving results")
        outname = os.path.split(filename)[1]
        outname = os.path.splitext(outname)[0]
        outname = "%s.sc" % outname
        fh = open(outname, 'w')
        fh.write("################################\n")
        fh.write("#                              #\n")
        fh.write("# Settings:                    #\n")
        fh.write("#  Method: point source        #\n")
        fh.write(f"#  Ref Ant: {args.reference:4d}               #\n")
        fh.write(f"#  Lower: {args.lower/1e6:5.1f} MHz            #\n")
        fh.write(f"#  Upper: {args.upper/1e6:5.1f} MHz            #\n")
        fh.write(f"#  Min (u,v): {args.min_uv_dist:4.1f} lambda      #\n")
        fh.write(f"#  Sun Factor: {args.sun_factor:6.0f}          #\n")
        fh.write(f"#  Max Iters: {args.max_iterations:3d}              #\n")
        fh.write(f"#  Delay Cutoff: {args.delay_cutoff:4.2f} ns       #\n")
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
        for i in range(delaysXX.size):
            fh.write("%3i  %.6g  %.6g  %.6g  %.6g\n" % (idi.stands[i], 1.0, delaysXX[i], 1.0, delaysYY[i]))
        fh.close()

        # Build up the images for each polarization
        if args.plot:
            print("    Gridding")
            toWork = numpy.where((freq>=80e6) & (freq<=82e6))[0]
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
                    #out = numpy.rot90(out)
                    #out = numpy.rot90(out)
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
                graticle(ax, lo.sidereal_time(), lo.lat)
                
            plt.show()
            
    print("...Done")


if __name__ == "__main__":
    numpy.seterr(all='ignore')
    
    parser = argparse.ArgumentParser(
        description="self-calibrate a TBF FITS IDI file",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
    parser.add_argument('filename', type=str, 
                        help='filename to calibrate')
    parser.add_argument('-r', '--reference', type=aph.positive_int, default=173,
                        help='reference stand to use')
    parser.add_argument('-l', '--lower', type=aph.positive_float, default=35.0,
                        help='lowest frequency to consider in MHz')
    parser.add_argument('-u', '--upper', type=aph.positive_float, default=85.0,
                        help='highest frequency to consider in MHz')
    parser.add_argument('-m', '--min-uv-dist', type=aph.positive_or_zero_float, default=14.0,
                        help='minimum baseline (u,v) length to use in wavelengths')
    parser.add_argument('-s', '--sun-factor', type=aph.positive_or_zero_float, default=1.0,
                        help="scale factor for the Sun's flux - useful for when it is flaring")
    parser.add_argument('-i', '--max-iterations', type=aph.positive_int, default=60,
                        help="maximum number of self-cal iterations")
    parser.add_argument('-d', '--delay-cutoff', type=aph.positive_float, default=0.2,
                        help="delay cutoff in ns for the self-cal convergence threshold")
    parser.add_argument('-p', '--plot', action='store_true',
                        help='plot the results at the end')
    args = parser.parse_args()
    args.lower *= 1e6
    args.upper *= 1e6
    main(args)
