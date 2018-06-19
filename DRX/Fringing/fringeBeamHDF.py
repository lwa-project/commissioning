#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Frank Schinzel's script to fringe special DRX files that have a beam X pol. 
and a dipole on Y pol.  The visibilities are written to an HDF file.

$Rev$
$LastChangedBy$
$LastChangedDate$
"""

import os
import sys
import numpy
import getopt
from datetime import datetime

from lsl.statistics import robust

from lsl.reader import drx
from lsl.common.dp import fS
from lsl.common import stations
from lsl.astro import unix_to_utcjd, DJD_OFFSET
from lsl.correlator import fx as fxc
from lsl.common.progress import ProgressBar

from matplotlib import pyplot as plt

import h5py


def usage(exitCode=None):
    print """fringeBeamHDF.py - Take a DRX file with the beam on X pol. and a dipole
on the Y pol. and cross correlate it.  The results are written out to an 
HDF5 file rather than a collection of .npz files.

Usage: fringeBeamHDF.py [OPTION] <dipole_ID_Y> <DRX_file>

Options:
-h, --help                  Display this help information
-l, --fft-length            Set FFT length (default = 512)
-t, --avg-time              Window to average visibilities in time (seconds; 
                            default = 4 s)
-d, --duration              Number of seconds to calculate the waterfall for 
                            (default = 10)
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
    config['duration'] = 0.0
    
    # Read in and process the command line flags
    try:
        opts, args = getopt.getopt(args, "hl:t:s:d:", ["help", "fft-length=", "avg-time=", "skip=", "duration="])
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
        elif opt in ('-d', '--duration'):
            config['duration'] = float(value)
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
    rawAntennas = site.getAntennas()
    
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
        nFramesFile = os.path.getsize(filename) / drx.FrameSize
        #junkFrame = drx.readFrame(fh)
        #fh.seek(0)
        while True:
            try:
                junkFrame = drx.readFrame(fh)
                try:
                    srate = junkFrame.getSampleRate()
                    t0 = junkFrame.getTime()
                    break
                except ZeroDivisionError:
                    pass
            except errors.syncError:
                fh.seek(-drx.FrameSize+1, 1)
                    
        fh.seek(-drx.FrameSize, 1)
        
        beam, tune, pol = junkFrame.parseID()
        srate = junkFrame.getSampleRate()
        
        tunepols = drx.getFramesPerObs(fh)
        tunepols = tunepols[0] + tunepols[1] + tunepols[2] + tunepols[3]
        beampols = tunepols
        
        # Offset in frames for beampols beam/tuning/pol. sets
        offset = int(config['offset'] * srate / 4096 * beampols)
        offset = int(1.0 * offset / beampols) * beampols
        fh.seek(offset*drx.FrameSize, 1)
        
        # Iterate on the offsets until we reach the right point in the file.  This
        # is needed to deal with files that start with only one tuning and/or a 
        # different sample rate.  
        while True:
            ## Figure out where in the file we are and what the current tuning/sample 
            ## rate is
            junkFrame = drx.readFrame(fh)
            srate = junkFrame.getSampleRate()
            t1 = junkFrame.getTime()
            tunepols = drx.getFramesPerObs(fh)
            tunepol = tunepols[0] + tunepols[1] + tunepols[2] + tunepols[3]
            beampols = tunepol
            fh.seek(-drx.FrameSize, 1)
            
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
            fh.seek(cOffset*drx.FrameSize, 1)
            
        # Update the offset actually used
        config['offset'] = t1 - t0
        offset = int(round(config['offset'] * srate / 4096 * beampols))
        offset = int(1.0 * offset / beampols) * beampols
        
        tnom = junkFrame.header.timeOffset
        tStart = junkFrame.getTime()
        
        # Get the DRX frequencies
        cFreq1 = 0.0
        cFreq2 = 0.0
        for i in xrange(4):
            junkFrame = drx.readFrame(fh)
            b,t,p = junkFrame.parseID()
            if p == 0 and t == 1:
                cFreq1 = junkFrame.getCentralFreq()
            elif p == 0 and t == 2:
                cFreq2 = junkFrame.getCentralFreq()
            else:
                pass
        fh.seek(-4*drx.FrameSize, 1)
        
        # Align the files as close as possible by the time tags and then make sure that
        # the first frame processed is from tuning 1, pol 0.
        junkFrame = drx.readFrame(fh)
        beam, tune, pol = junkFrame.parseID()
        pair = 2*(tune-1) + pol
        j = 0
        while pair != 0:
            junkFrame = drx.readFrame(fh)
            beam, tune, pol = junkFrame.parseID()
            pair = 2*(tune-1) + pol
            j += 1
        fh.seek(-drx.FrameSize, 1)
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
        print "  Duration of File: %f" % tFile
        print "  Offset: %f s" % offset
        
        if config['duration']!=0:
            nChunks = int(round(config['duration'] / tInt))
        else:
            nChunks = int(tFile/tInt)
            
        print "Processing: %i integrations" % nChunks
        
        # Here we start the HDF5 file
        outname = os.path.split(filename)[1]
        outname = os.path.splitext(outname)[0]
        outname = "%s.hdf5" % outname
        outfile = h5py.File(outname)
        group1 = outfile.create_group("Time")
        group2 = outfile.create_group("Frequencies")
        group3 = outfile.create_group("Visibilities")
        out = raw_input("Target Name: ")
        outfile.attrs["OBJECT"] = out
        out = raw_input("Polarization (X/Y): ")
        outfile.attrs["POLARIZATION"] = out
        dset1 = group1.create_dataset("Timesteps", (nChunks,3), numpy.float64, maxshape=(nChunks,3))
        dset2 = group2.create_dataset("Tuning1",(LFFT-1 if float(fxc.__version__) < 0.8 else LFFT, ), '<f8', maxshape=(LFFT-1 if float(fxc.__version__) < 0.8 else LFFT, ) )
        dset3 = group2.create_dataset("Tuning2",(LFFT-1 if float(fxc.__version__) < 0.8 else LFFT, ), '<f8', maxshape=(LFFT-1 if float(fxc.__version__) < 0.8 else LFFT, ) )
        dset4 = group3.create_dataset("Tuning1", (nChunks, 3, LFFT-1 if float(fxc.__version__) < 0.8 else LFFT), numpy.complex64,maxshape=(nChunks, 3, LFFT-1 if float(fxc.__version__) < 0.8 else LFFT))
        dset5 = group3.create_dataset("Tuning2", (nChunks, 3, LFFT-1 if float(fxc.__version__) < 0.8 else LFFT), numpy.complex64,maxshape=(nChunks, 3, LFFT-1 if float(fxc.__version__) < 0.8 else LFFT))
        
        pb = ProgressBar(max=nChunks)
        tsec = numpy.zeros(1, dtype=numpy.float64)
        for i in xrange(nChunks):
            junkFrame = drx.readFrame(fh)
            tStart = junkFrame.getTime()
            fh.seek(-drx.FrameSize, 1)
            
            count1 = [0,0]
            data1 = numpy.zeros((2, 4096*nFrames), dtype=numpy.complex64)
            count2 = [0,0]
            data2 = numpy.zeros((2, 4096*nFrames), dtype=numpy.complex64)
            for j in xrange(nFrames):
                for k in xrange(4):
                    cFrame = drx.readFrame(fh)
                    beam, tune, pol = cFrame.parseID()
                    pair = 2*(tune-1) + pol
                    
                    if tune == 1:
                        data1[pol, count1[pol]*4096:(count1[pol]+1)*4096] = cFrame.data.iq
                        count1[pol] += 1
                    else:
                        data2[pol, count2[pol]*4096:(count2[pol]+1)*4096] = cFrame.data.iq
                        count2[pol] += 1
                        
            # Correlate
            blList1, freq1, vis1 = fxc.FXMaster(data1, antennas, LFFT=LFFT, Overlap=1, IncludeAuto=True, verbose=False, SampleRate=srate, CentralFreq=cFreq1, Pol='XX', ReturnBaselines=True, GainCorrect=False, ClipLevel=0)
            
            blList2, freq2, vis2 = fxc.FXMaster(data2, antennas, LFFT=LFFT, Overlap=1, IncludeAuto=True, verbose=False, SampleRate=srate, CentralFreq=cFreq2, Pol='XX', ReturnBaselines=True, GainCorrect=False, ClipLevel=0)
            
            if i == 0:
                tsec = tInt/2
                outfile.attrs["STANDS"] = numpy.array([stand1, stand2])
                outfile.attrs["SRATE"] = srate
                date = datetime.fromtimestamp(tStart).date()
                outfile.attrs["DATE"] = str(date)
                dset2.write_direct(freq1)
                dset3.write_direct(freq2)
            else:
                tsec += tInt
                
            temp = numpy.zeros(3, dtype=numpy.float64)
            temp[0] = tStart
            temp[1] = tInt
            temp[2] = tsec
            dset1.write_direct(temp, dest_sel=numpy.s_[i])
            dset4.write_direct(vis1, dest_sel=numpy.s_[i])
            dset5.write_direct(vis2, dest_sel=numpy.s_[i])
            del data1
            del data2
            
            pb.inc(amount=1)
            sys.stdout.write(pb.show()+'\r')
            sys.stdout.flush()
            
        sys.stdout.write(pb.show()+'\r')
        sys.stdout.write('\n')
        sys.stdout.flush()
        outfile.close()
        
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
    
