#!/usr/bin/env python3

"""
Script to fringe special DRX files that have a beam X pol. and a dipole 
on Y pol.  The visibilites are written to a NPZ file.
"""

import os
import sys
import numpy
import argparse

from lsl.statistics import robust

from lsl.reader import drx, errors, buffer
from lsl.common import stations
from lsl.astro import unix_to_utcjd, DJD_OFFSET
from lsl.correlator import fx as fxc
from lsl.common.progress import ProgressBarPlus
from lsl.misc import parser as aph

from matplotlib import pyplot as plt


def main(args):
    LFFT = args.fft_length
    
    stand1 = 0
    stand2 = int(args.dipole_id_y)
    filenames = args.filename
    
    # Build up the station
    if args.lwasv:
        site = stations.lwasv
    else:
        site = stations.lwa1
        
    # Get the antennas we need (and a fake one for the beam)
    rawAntennas = site.antennas
    
    antennas = []
    
    dipole = None
    xyz = numpy.zeros((len(rawAntennas),3))
    i = 0
    for ant in rawAntennas:
        if ant.stand.id == stand2 and ant.pol == 0:
            dipole = ant
        xyz[i,0] = ant.stand.x
        xyz[i,1] = ant.stand.y
        xyz[i,2] = ant.stand.z
        i += 1
    arrayX = xyz[:,0].mean()
    arrayY = xyz[:,1].mean()
    arrayZ = xyz[:,2].mean()
    
    ## Fake one down here...
    beamStand   = stations.Stand(0, arrayX, arrayY, arrayZ)
    beamFEE     = stations.FEE('Beam', 0, gain1=0, gain2=0, status=3)
    beamCable   = stations.Cable('Beam', 0, vf=1.0)
    beamAntenna = stations.Antenna(0, stand=beamStand, pol=0, theta=0, phi=0, status=3)
    beamAntenna.fee = beamFEE
    beamAntenna.feePort = 1
    beamAntenna.cable = beamCable
    
    antennas.append(beamAntenna)
    
    ## Dipole down here...
    ### NOTE
    ### Here we zero out the cable length for the dipole since the delay 
    ### setup that is used for these observations already takes the 
    ### cable/geometric delays into account.  We shouldn't need anything 
    ### else to get good fringes.
    dipole.cable.length = 0
    antennas.append(dipole)
    
    # Loop over the input files...
    for filename in filenames:
        fh = open(filename, "rb")
        nFramesFile = os.path.getsize(filename) // drx.FRAME_SIZE
        #junkFrame = drx.read_frame(fh)
        #fh.seek(0)
        while True:
            try:
                junkFrame = drx.read_frame(fh)
                try:
                    srate = junkFrame.sample_rate
                    t0 =  junkFrame.time
                    break
                except ZeroDivisionError:
                    pass
            except errors.SyncError:
                fh.seek(-drx.FRAME_SIZE+1, 1)
                
        fh.seek(-drx.FRAME_SIZE, 1)
        
        beam, tune, pol = junkFrame.id
        srate = junkFrame.sample_rate
        
        tunepols = drx.get_frames_per_obs(fh)
        tunepols = tunepols[0] + tunepols[1] + tunepols[2] + tunepols[3]
        beampols = tunepols
        
        # Offset in frames for beampols beam/tuning/pol. sets
        offset = int(args.skip * srate / 4096 * beampols)
        offset = int(1.0 * offset / beampols) * beampols
        fh.seek(offset*drx.FRAME_SIZE, 1)
        
        # Iterate on the offsets until we reach the right point in the file.  This
        # is needed to deal with files that start with only one tuning and/or a 
        # different sample rate.  
        while True:
            ## Figure out where in the file we are and what the current tuning/sample 
            ## rate is
            junkFrame = drx.read_frame(fh)
            srate = junkFrame.sample_rate
            t1 = junkFrame.time
            tunepols = drx.get_frames_per_obs(fh)
            tunepol = tunepols[0] + tunepols[1] + tunepols[2] + tunepols[3]
            beampols = tunepol
            fh.seek(-drx.FRAME_SIZE, 1)
            
            ## See how far off the current frame is from the target
            tDiff = t1 - (t0 + args.skip)
            
            ## Half that to come up with a new seek parameter
            tCorr = -tDiff / 8.0
            cOffset = int(tCorr * srate / 4096 * beampols)
            cOffset = int(1.0 * cOffset / beampols) * beampols
            offset += cOffset
            
            ## If the offset is zero, we are done.  Otherwise, apply the offset
            ## and check the location in the file again/
            if cOffset is 0:
                    break
            fh.seek(cOffset*drx.FRAME_SIZE, 1)
            
        # Update the offset actually used
        args.skip = t1 - t0
        offset = int(round(args.skip * srate / 4096 * beampols))
        offset = int(1.0 * offset / beampols) * beampols
        
        tnom = junkFrame.header.time_offset
        tStart = junkFrame.time
        
        # Get the DRX frequencies
        cFreq1 = 0.0
        cFreq2 = 0.0
        for i in range(32):
            junkFrame = drx.read_frame(fh)
            b,t,p = junkFrame.id
            if p == 0 and t == 1:
                cFreq1 = junkFrame.central_freq
            elif p == 0 and t == 2:
                cFreq2 = junkFrame.central_freq
            else:
                pass
        fh.seek(-32*drx.FRAME_SIZE, 1)
        
        # Align the files as close as possible by the time tags and then make sure that
        # the first frame processed is from tuning 1, pol 0.
        junkFrame = drx.read_frame(fh)
        beam, tune, pol = junkFrame.id
        pair = 2*(tune-1) + pol
        j = 0
        while pair != 0:
            junkFrame = drx.read_frame(fh)
            beam, tune, pol = junkFrame.id
            pair = 2*(tune-1) + pol
            j += 1
        fh.seek(-drx.FRAME_SIZE, 1)
        print("Shifted beam %i data by %i frames (%.4f s)" % (beam, j, j*4096/srate/4))
        
        # Set integration time
        tInt = args.avg_time
        nFrames = int(round(tInt*srate/4096))
        tInt = nFrames*4096/srate
        
        # Read in some data
        tFile = nFramesFile / 4 * 4096 / srate
        
        # Report
        print("Filename: %s" % filename)
        print("  Sample Rate: %i Hz" % srate)
        print("  Tuning 1: %.1f Hz" % cFreq1)
        print("  Tuning 2: %.1f Hz" % cFreq2)
        print("  ===")
        print("  Integration Time: %.3f s" % tInt)
        print("  Integrations in File: %i" % int(tFile/tInt))

        nChunks = int(tFile/tInt)
        drxBuffer = buffer.DRXFrameBuffer(beams=[beam,], tunes=[1,2], pols=[0,1])
         
        data = numpy.zeros((2,2,4096*nFrames), dtype=numpy.complex64)

        pb = ProgressBarPlus(max=nChunks)
        for i in range(nChunks):
            j = 0
            while j < nFrames:
                for k in range(4):
                    try:
                        cFrame = drx.read_frame(fh)
                        drxBuffer.append( cFrame )
                    except errors.SyncError:
                        pass
    
                cFrames = drxBuffer.get()
                if cFrames is None:
                    continue

                for cFrame in cFrames:
                    if j == 0:
                        tStart = cFrame.time
                    beam, tune, pol = cFrame.id
                    pair = 2*(tune-1) + pol

                    if tune == 1:
                        data[0, pol, j*4096:(j+1)*4096] = cFrame.payload.data
                    else:
                        data[1, pol, j*4096:(j+1)*4096] = cFrame.payload.data

                j += 1
 
            # Correlate
            blList1, freq1, vis1 = fxc.FXMaster(data[0,:,:], antennas, LFFT=LFFT, overlap=1, include_auto=True, verbose=False, sample_rate=srate, central_freq=cFreq1, pol='XX', return_baselines=True, gain_correct=False, clip_level=0)
        
            blList2, freq2, vis2 = fxc.FXMaster(data[1,:,:], antennas, LFFT=LFFT, overlap=1, include_auto=True, verbose=False, sample_rate=srate, central_freq=cFreq2, pol='XX', return_baselines=True, gain_correct=False, clip_level=0)
            
            if nChunks != 1:
                outfile = os.path.split(filename)[1]
                outfile = os.path.splitext(outfile)[0]
                outfile = "%s-vis-%04i.npz" % (outfile, i+1)
            else:
                outfile = os.path.split(filename)[1]
                outfile = os.path.splitext(outfile)[0]
                outfile = "%s-vis.npz" % outfile
            numpy.savez(outfile, srate=srate, freq1=freq1, vis1=vis1, freq2=freq2, vis2=vis2, tStart=tStart.unix, tInt=tInt, stands=numpy.array([stand1, stand2]))
             
            pb.inc(amount=1)
            sys.stdout.write(pb.show()+'\r')
            sys.stdout.flush()
            
        sys.stdout.write(pb.show()+'\r')
        sys.stdout.write('\n')
        sys.stdout.flush()
        
        # Plot
        fig = plt.figure()
        i = 0
        for bl, vi in zip(blList1, vis1):
            ax = fig.add_subplot(4, 3, i+1)
            ax.plot(freq1/1e6, numpy.unwrap(numpy.angle(vi)))
            ax.set_title('Stand %i - Stand %i' % (bl[0].stand.id, bl[1].stand.id))
            ax = fig.add_subplot(4, 3, i+4)
            ax.plot(freq1/1e6, numpy.abs(vi))
            i += 1
            
            coeff = numpy.polyfit(freq1, numpy.unwrap(numpy.angle(vi)), 1)
            #print(coeff[0]/2/numpy.pi*1e9, coeff[1]*180/numpy.pi)
            
        i = 6
        for bl, vi in zip(blList2, vis2):
            ax = fig.add_subplot(4, 3, i+1)
            ax.plot(freq2/1e6, numpy.unwrap(numpy.angle(vi)))
            ax.set_title('Stand %i - Stand %i' % (bl[0].stand.id, bl[1].stand.id))
            ax = fig.add_subplot(4, 3, i+4)
            ax.plot(freq2/1e6, numpy.abs(vi))
            i += 1
            
            coeff = numpy.polyfit(freq2, numpy.unwrap(numpy.angle(vi)), 1)
            #print(coeff[0]/2/numpy.pi*1e9, coeff[1]*180/numpy.pi)
            
        #plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="take a DRX file with the beam on X pol. and a dipole on the Y pol. and cross correlate it",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
    parser.add_argument('dipole_id_y', type=aph.positive_int,
                        help='dipole number')
    parser.add_argument('filename', type=str, nargs='+',
                        help='filename to process')
    parser.add_argument('-v', '--lwasv', action='store_true',
                        help='Station is LWA-SV')
    parser.add_argument('-l', '--fft-length', type=aph.positive_int, default=512,
                        help='FFT transform size')
    parser.add_argument('-t', '--avg-time', type=aph.positive_float, default=4.0,
                        help='window to average visibilities in seconds')
    parser.add_argument('-s', '--skip', type=aph.positive_or_zero_float, default=0.0,
                        help='skip the specified number of seconds at the beginning of the file')
    args = parser.parse_args()
    main(args)
    
