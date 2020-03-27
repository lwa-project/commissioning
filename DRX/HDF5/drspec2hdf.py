#!/usr/bin/env python

"""
Small script to read in a DR spectrometer binary data file and create a HDF5 in 
the image of hdfWaterfall.py that can be plotted with plotHDF.py
"""

# Python3 compatiability
from __future__ import print_function, division
import sys
if sys.version_info > (3,):
    xrange = range
    raw_input = input
    
import os
import sys
import h5py
import numpy
import argparse
from datetime import datetime

from lsl.reader import drx, drspec, errors
from lsl.common import progress
from lsl.common import mcs, sdf, metabundle
from lsl.misc import parser as aph
try:
    from lsl.common import sdfADP, metabundleADP
    adpReady = True
except ImportError:
    adpReady = False

import data as hdfData


def main(args):
    # Set the site
    site = None
    if args.lwa1:
        site = 'lwa1'
    elif args.lwasv:
        site = 'lwasv'
        
    # Open the file and file good data (not raw DRX data)
    fh = open(args.filename, 'rb')

    try:
        for i in xrange(5):
            junkFrame = drx.read_frame(fh)
        raise RuntimeError("ERROR: '%s' appears to be a raw DRX file, not a DR spectrometer file" % args.filename)
    except errors.SyncError:
        fh.seek(0)
        
    # Interrogate the file to figure out what frames sizes to expect, now many 
    # frames there are, and what the transform length is
    FRAME_SIZE = drspec.get_frame_size(fh)
    nFrames = os.path.getsize(args.filename) // FRAME_SIZE
    nChunks = nFrames
    LFFT = drspec.get_transform_size(fh)

    # Read in the first frame to figure out the DP information
    junkFrame = drspec.read_frame(fh)
    fh.seek(-FRAME_SIZE, 1)
    srate = junkFrame.sample_rate
    t0 = junkFrame.get_time()
    tInt = junkFrame.header.nInts*LFFT/srate
    
    # Offset in frames for beampols beam/tuning/pol. sets
    offset = int(round(args.skip / tInt))
    fh.seek(offset*FRAME_SIZE, 1)
    
    # Iterate on the offsets until we reach the right point in the file.  This
    # is needed to deal with files that start with only one tuning and/or a 
    # different sample rate.  
    while True:
        ## Figure out where in the file we are and what the current tuning/sample 
        ## rate is
        junkFrame = drspec.read_frame(fh)
        srate = junkFrame.sample_rate
        t1 = junkFrame.get_time()
        tInt = junkFrame.header.nInts*LFFT/srate
        fh.seek(-FRAME_SIZE, 1)
        
        ## See how far off the current frame is from the target
        tDiff = t1 - (t0 + args.skip)
        
        ## Half that to come up with a new seek parameter
        tCorr = -tDiff / 2.0
        cOffset = int(round(tCorr / tInt))
        offset += cOffset
        
        ## If the offset is zero, we are done.  Otherwise, apply the offset
        ## and check the location in the file again/
        if cOffset is 0:
            break
        fh.seek(cOffset*FRAME_SIZE, 1)
        
    # Update the offset actually used
    args.skip = t1 - t0
    nChunks = (os.path.getsize(args.filename) - fh.tell()) // FRAME_SIZE
    
    # Update the file contents
    beam = junkFrame.id
    central_freq1 = junkFrame.get_central_freq(1)
    central_freq2 = junkFrame.get_central_freq(2)
    srate = junkFrame.sample_rate
    data_products = junkFrame.data_products
    t0 = junkFrame.get_time()
    tInt = junkFrame.header.nInts*LFFT/srate
    beginDate = datetime.utcfromtimestamp(junkFrame.get_time())
        
    # Report
    print("Filename: %s" % args.filename)
    if args.metadata is not None:
        print("Metadata: %s" % args.metadata)
    elif args.sdf is not None:
        print("SDF: %s" % args.sdf)
    print("Date of First Frame: %s" % beginDate)
    print("Beam: %i" % beam)
    print("Sample Rate: %i Hz" % srate)
    print("Tuning Frequency: %.3f Hz (1); %.3f Hz (2)" % (central_freq1, central_freq2))
    print("Data Products: %s" % ','.join(data_products))
    print("Frames: %i (%.3f s)" % (nFrames, nFrames*tInt))
    print("---")
    print("Offset: %.3f s (%i frames)" % (args.skip, offset))
    print("Transform Length: %i" % LFFT)
    print("Integration: %.3f s" % tInt)
    
    # Setup the output file
    outname = os.path.split(args.filename)[1]
    outname = os.path.splitext(outname)[0]
    outname = '%s-waterfall.hdf5' % outname
    
    if os.path.exists(outname):
        if not args.force:
            yn = raw_input("WARNING: '%s' exists, overwrite? [Y/n] " % outname)
        else:
            yn = 'y'
            
        if yn not in ('n', 'N'):
            os.unlink(outname)
        else:
            raise RuntimeError("Output file '%s' already exists" % outname)
            
    f = hdfData.createNewFile(outname)
    obsList = {}
    if args.metadata is not None:
        try:
            project = metabundle.get_sdf(args.metadata)
        except Exception as e:
            if adpReady:
                project = metabundleADP.get_sdf(args.metadata)
            else:
                raise e
                
        sdfBeam  = project.sessions[0].drxBeam
        spcSetup = project.sessions[0].spcSetup
        if sdfBeam != beam:
            raise RuntimeError("Metadata is for beam #%i, but data is from beam #%i" % (sdfBeam, beam))
            
        for i,obs in enumerate(project.sessions[0].observations):
            sdfStart = mcs.mjdmpm_to_datetime(obs.mjd, obs.mpm)
            sdfStop  = mcs.mjdmpm_to_datetime(obs.mjd, obs.mpm + obs.dur)
            obsChunks = int(numpy.ceil(obs.dur/1000.0 * drx.FILTER_CODES[obs.filter] / (spcSetup[0]*spcSetup[1])))
            
            obsList[i+1] = (sdfStart, sdfStop, obsChunks)
            
        hdfData.fillFromMetabundle(f, args.metadata)
        
    elif args.sdf is not None:
        try:
            project = sdf.parse_sdf(args.sdf)
        except Exception as e:
            if adpReady:
                project = sdfADP.parse_sdf(args.sdf)
            else:
                raise e
                
        sdfBeam  = project.sessions[0].drxBeam
        spcSetup = project.sessions[0].spcSetup
        if sdfBeam != beam:
            raise RuntimeError("Metadata is for beam #%i, but data is from beam #%i" % (sdfBeam, beam))
            
        for i,obs in enumerate(project.sessions[0].observations):
            sdfStart = mcs.mjdmpm_to_datetime(obs.mjd, obs.mpm)
            sdfStop  = mcs.mjdmpm_to_datetime(obs.mjd, obs.mpm + obs.dur)
            obsChunks = int(numpy.ceil(obs.dur/1000.0 * drx.FILTER_CODES[obs.filter] / (spcSetup[0]*spcSetup[1])))
            
            obsList[i+1] = (sdfStart, sdfStop, obsChunks)
            
        hdfData.fillFromSDF(f, args.sdf, station=site)
        
    else:
        obsList[1] = (beginDate, datetime(2222,12,31,23,59,59), nChunks)
        
        hdfData.fillMinimum(f, 1, beam, srate, station=site)
        
    data_products = junkFrame.data_products
    for o in sorted(obsList.keys()):
        for t in (1,2):
            hdfData.createDataSets(f, o, t, numpy.arange(LFFT, dtype=numpy.float64), obsList[o][2], data_products)
            
    f.attrs['FileGenerator'] = 'drspec2hdf.py'
    f.attrs['InputData'] = os.path.basename(args.filename)
    
    # Create the various HDF group holders
    ds = {}
    for o in sorted(obsList.keys()):
        obs = hdfData.getObservationSet(f, o)
        
        ds['obs%i' % o] = obs
        ds['obs%i-time' % o] = obs.create_dataset('time', (obsList[o][2],), 'f8')
        
        for t in (1,2):
            ds['obs%i-freq%i' % (o, t)] = hdfData.get_data_set(f, o, t, 'freq')
            for p in data_products:
                ds["obs%i-%s%i" % (o, p, t)] = hdfData.get_data_set(f, o, t, p)
            ds['obs%i-Saturation%i' % (o, t)] = hdfData.get_data_set(f, o, t, 'Saturation')
            
    # Loop over DR spectrometer frames to fill in the HDF5 file
    pbar = progress.ProgressBar(max=nChunks)
    o = 1
    j = 0
    
    firstPass = True
    for i in xrange(nChunks):
        frame = drspec.read_frame(fh)
        
        cTime = datetime.utcfromtimestamp(frame.get_time())
        if cTime > obsList[o][1]:
            # Increment to the next observation
            o += 1
            
            # If we have reached the end, exit...
            try:
                obsList[o]
                
                firstPass = True
            except KeyError:
                sys.stdout.write('%s\r' % (' '*pbar.span))
                sys.stdout.flush()
                print("End of observing block according to SDF, exiting")
                break
                
        if cTime < obsList[o][0]:
            # Skip over data that occurs before the start of the observation
            continue
            
        try:
            if frame.get_time() > oTime + 1.001*tInt:
                print('Warning: Time tag error at frame %i; %.3f > %.3f + %.3f' % (i, frame.get_time(), oTime, tInt))
        except NameError:
            pass
        oTime = frame.get_time()
        
        if firstPass:
            # Otherwise, continue on...
            central_freq1 = frame.get_central_freq(1)
            central_freq2 = frame.get_central_freq(2)
            srate = frame.sample_rate
            tInt  = frame.header.nInts*LFFT/srate
            
            freq = numpy.fft.fftshift( numpy.fft.fftfreq(LFFT, d=1.0/srate) )
            freq = freq.astype(numpy.float64)
            
            sys.stdout.write('%s\r' % (' '*pbar.span))
            sys.stdout.flush()
            print("Switching to Obs. #%i" % o)
            print("-> Tunings: %.1f Hz, %.1f Hz" % (central_freq1, central_freq2))
            print("-> Sample Rate: %.1f Hz" % srate)
            print("-> Integration Time: %.3f s" % tInt)
            sys.stdout.write(pbar.show()+'\r')
            sys.stdout.flush()
            
            j = 0
            ds['obs%i-freq1' % o][:] = freq + central_freq1
            ds['obs%i-freq2' % o][:] = freq + central_freq2
            
            obs = ds['obs%i' % o]
            obs.attrs['tInt'] = tInt
            obs.attrs['tInt_Units'] = 's'
            obs.attrs['LFFT'] = LFFT
            obs.attrs['nchan'] = LFFT
            obs.attrs['RBW'] = freq[1]-freq[0]
            obs.attrs['RBW_Units'] = 'Hz'
            
            firstPass = False
            
        # Load the data from the spectrometer frame into the HDF5 group
        ds['obs%i-time' % o][j] = frame.get_time()
        
        ds['obs%i-Saturation1' % o][j,:] = frame.data.saturations[0:2]
        ds['obs%i-Saturation2' % o][j,:] = frame.data.saturations[2:4]
        
        for t in (1,2):
            for p in data_products:
                ds['obs%i-%s%i' % (o, p, t)][j,:] = getattr(frame.data, "%s%i" % (p, t-1), None)
        j += 1
        
        # Update the progress bar
        pbar.inc()
        if i % 10 == 0:
            sys.stdout.write(pbar.show()+'\r')
            sys.stdout.flush()
            
    sys.stdout.write(pbar.show()+'\n')
    sys.stdout.flush()
    
    # Done
    fh.close()

    # Save the output to a HDF5 file
    f.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='convert a DR spectrometer file to a HDF5 similar to what hdfWaterfall.py generates', 
        epilog='NOTE:  Both the -m/--metadata and -d/--sdf options provide the same additional observation information to drspec2hdf.py so only one needs to be provided.', 
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
    parser.add_argument('filename', type=str, 
                        help='filename to convert')
    parser.add_argument('-s', '--skip', type=aph.positive_or_zero_float, default=0.0, 
                        help='skip foward the specified number of seconds into the file')
    parser.add_argument('-m', '--metadata', type=str, 
                        help='metadata tarball for additional information')
    parser.add_argument('-d', '--sdf', type=str, 
                        help='SDF for additional information')
    sgroup = parser.add_mutually_exclusive_group(required=False)
    sgroup.add_argument('-1', '--lwa1', action='store_true', default=True, 
                        help='data is from LWA1; needed for -d/--sdf or when no metadata is provided')
    sgroup.add_argument('-v', '--lwasv', action='store_true', 
                        help='data is from LWA-SV; needed for -d/--sdf or when no metadata is provided')
    parser.add_argument('-f', '--force', action='store_true', 
                        help='force overwritting of existing HDF5 files')
    args = parser.parse_args()
    main(args)
    