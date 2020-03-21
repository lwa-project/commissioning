#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Script to fringe special DRX files that have a beam X pol. and a dipole 
on Y pol.  The visibilites are written to a NPZ file.
"""

import os
import sys
import numpy
import getopt

from lsl.statistics import robust

from lsl.reader import drx, errors
from lsl.common import stations
from lsl.astro import unix_to_utcjd, DJD_OFFSET
from lsl.correlator import fx as fxc
from lsl.common.progress import ProgressBar

from matplotlib import pyplot as plt


def usage(exitCode=None):
    print """fringeBeam.py - Take a DRX file with the beam on X pol. and a dipole
on the Y pol. and cross correlate it.

Usage: fringeBeam.py [OPTION] <dipole_ID_Y> <DRX_file>

Options:
-h, --help                  Display this help information
-l, --fft-length            Set FFT length (default = 512)
-t, --avg-time              Window to average visibilities in time (seconds; 
                            default = 4 s)
-s, --skip                  Skip the specified number of seconds at the beginning
                            of the file (default = 0)
"""
    
    if exitCode is not None:
        sys.exit(exitCode)
    else:
        return True


def parseOptions(args):
    config = {}
    config['avgTime'] = 4.0
    config['LFFT'] = 512
    config['offset'] = 0.0
    
    # Read in and process the command line flags
    try:
        opts, args = getopt.getopt(args, "hl:t:s:", ["help", "fft-length=", "avg-time=", "skip="])
    except getopt.GetoptError, err:
        # Print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        usage(exitCode=2)
        
    # Work through opts
    for opt, value in opts:
        if opt in ('-h', '--help'):
            usage(exitCode=0)
        elif opt in ('-l', '--fft-length'):
            config['LFFT'] = int(value)
        elif opt in ('-t', '--avg-time'):
            config['avgTime'] = float(value)
        elif opt in ('-s', '--skip'):
            config['offset'] = float(value)
        else:
            assert False
            
    # Add in arguments
    config['args'] = args
    
    # Return configuration
    return config


def main(args):
    config = parseOptions(args)
    
    LFFT = config['LFFT']
    
    stand1 = 0
    stand2 = int(config['args'][0])
    filenames = config['args'][1:]
    
    # Build up the station
    site = stations.lwa1
    
    # Get the antennas we need (and a fake one for the beam)
    rawAntennas = site.antennas
    
    dipole = None
    for ant in rawAntennas:
        if ant.stand.id == stand2 and ant.pol == 0:
            dipole = ant
            
    antennas = []
    
    ## Fake one down here...
    beamStand   = stations.Stand(0, dipole.x, dipole.y, dipole.z)
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
        nFramesFile = os.path.getsize(filename) / drx.FRAME_SIZE
        #junkFrame = drx.read_frame(fh)
        #fh.seek(0)
        while True:
            try:
                junkFrame = drx.read_frame(fh)
                try:
                    srate = junkFrame.sample_rate
                    t0 = junkFrame.get_time()
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
        offset = int(config['offset'] * srate / 4096 * beampols)
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
            t1 = junkFrame.get_time()
            tunepols = drx.get_frames_per_obs(fh)
            tunepol = tunepols[0] + tunepols[1] + tunepols[2] + tunepols[3]
            beampols = tunepol
            fh.seek(-drx.FRAME_SIZE, 1)
            
            ## See how far off the current frame is from the target
            tDiff = t1 - (t0 + config['offset'])
            
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
        config['offset'] = t1 - t0
        offset = int(round(config['offset'] * srate / 4096 * beampols))
        offset = int(1.0 * offset / beampols) * beampols
        
        tnom = junkFrame.header.time_offset
        tStart = junkFrame.get_time()
        
        # Get the DRX frequencies
        cFreq1 = 0.0
        cFreq2 = 0.0
        for i in xrange(4):
            junkFrame = drx.read_frame(fh)
            b,t,p = junkFrame.id
            if p == 0 and t == 1:
                cFreq1 = junkFrame.central_freq
            elif p == 0 and t == 2:
                cFreq2 = junkFrame.central_freq
            else:
                pass
        fh.seek(-4*drx.FRAME_SIZE, 1)
        
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
        print "Shifted beam %i data by %i frames (%.4f s)" % (beam, j, j*4096/srate/4)
        
        # Set integration time
        tInt = config['avgTime']
        nFrames = int(round(tInt*srate/4096))
        tInt = nFrames*4096/srate
        
        # Read in some data
        tFile = nFramesFile / 4 * 4096 / srate
        
        # Report
        print "Filename: %s" % filename
        print "  Sample Rate: %i Hz" % srate
        print "  Tuning 1: %.1f Hz" % cFreq1
        print "  Tuning 2: %.1f Hz" % cFreq2
        print "  ==="
        print "  Integration Time: %.3f s" % tInt
        print "  Integrations in File: %i" % int(tFile/tInt)
        
        nChunks = int(tFile/tInt)
        pb = ProgressBar(max=nChunks)
        for i in xrange(nChunks):
            junkFrame = drx.read_frame(fh)
            tStart = junkFrame.get_time()
            fh.seek(-drx.FRAME_SIZE, 1)
            
            count1 = [0,0]
            data1 = numpy.zeros((2, 4096*nFrames), dtype=numpy.complex64)
            count2 = [0,0]
            data2 = numpy.zeros((2, 4096*nFrames), dtype=numpy.complex64)
            for j in xrange(nFrames):
                for k in xrange(4):
                    cFrame = drx.read_frame(fh)
                    beam, tune, pol = cFrame.id
                    pair = 2*(tune-1) + pol
                    
                    if tune == 1:
                        data1[pol, count1[pol]*4096:(count1[pol]+1)*4096] = cFrame.payload.data
                        count1[pol] += 1
                    else:
                        data2[pol, count2[pol]*4096:(count2[pol]+1)*4096] = cFrame.payload.data
                        count2[pol] += 1
                    
            # Correlate
            blList1, freq1, vis1 = fxc.FXMaster(data1, antennas, LFFT=LFFT, Overlap=1, include_auto=True, verbose=False, sample_rate=srate, central_freq=cFreq1, Pol='XX', return_baselines=True, gain_correct=False, clip_level=0)
        
            blList2, freq2, vis2 = fxc.FXMaster(data2, antennas, LFFT=LFFT, Overlap=1, include_auto=True, verbose=False, sample_rate=srate, central_freq=cFreq2, Pol='XX', return_baselines=True, gain_correct=False, clip_level=0)
            
            if nChunks != 1:
                outfile = os.path.split(filename)[1]
                outfile = os.path.splitext(outfile)[0]
                outfile = "%s-vis-%04i.npz" % (outfile, i+1)
            else:
                outfile = os.path.split(filename)[1]
                outfile = os.path.splitext(outfile)[0]
                outfile = "%s-vis.npz" % outfile
            numpy.savez(outfile, srate=srate, freq1=freq1, vis1=vis1, freq2=freq2, vis2=vis2, tStart=tStart, tInt=tInt, stands=numpy.array([stand1, stand2]))
            
            del data1
            del data2
            
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
            #print coeff[0]/2/numpy.pi*1e9, coeff[1]*180/numpy.pi
            
        i = 6
        for bl, vi in zip(blList2, vis2):
            ax = fig.add_subplot(4, 3, i+1)
            ax.plot(freq2/1e6, numpy.unwrap(numpy.angle(vi)))
            ax.set_title('Stand %i - Stand %i' % (bl[0].stand.id, bl[1].stand.id))
            ax = fig.add_subplot(4, 3, i+4)
            ax.plot(freq2/1e6, numpy.abs(vi))
            i += 1
            
            coeff = numpy.polyfit(freq2, numpy.unwrap(numpy.angle(vi)), 1)
            #print coeff[0]/2/numpy.pi*1e9, coeff[1]*180/numpy.pi
            
        #plt.show()


if __name__ == "__main__":
    main(sys.argv[1:])
    
