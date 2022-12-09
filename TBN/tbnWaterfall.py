#!/usr/bin/env python3

"""
Given a TBN file, plot the time averaged spectra for each stand output over some 
period.
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
import math
import numpy
import argparse
from datetime import datetime

from lsl.reader import tbn, errors
from lsl.reader.ldp import LWA1DataFile
import lsl.correlator.fx as fxc
from lsl.astro import unix_to_utcjd, DJD_OFFSET
from lsl.common import progress, stations
from lsl.common import mcs, sdf, metabundle
from lsl.common import sdfADP, metabundleADP
from lsl.misc import parser as aph

import matplotlib.pyplot as plt

import data as hdfData


def best_freq_units(freq):
    """Given a numpy array of frequencies in Hz, return a new array with the
    frequencies in the best units possible (kHz, MHz, etc.)."""
    
    # Figure out how large the data are
    scale = int(math.log10(freq.max()))
    if scale >= 9:
        divis = 1e9
        units = 'GHz'
    elif scale >= 6:
        divis = 1e6
        units = 'MHz'
    elif scale >= 3:
        divis = 1e3
        units = 'kHz'
    else:
        divis = 1
        units = 'Hz'
        
    # Convert the frequency
    newFreq = freq / divis
    
    # Return units and freq
    return (newFreq, units)


def process_data_to_linear(idf, antennas, tStart, duration, sample_rate, args, dataSets, obsID=1, clip1=0):
    """
    Process a chunk of data in a raw DRX file into linear polarization 
    products and add the contents to an HDF5 file.
    """
    
    # Length of the FFT
    LFFT = args.fft_length
    
    # Find the start of the observation
    print('Looking for #%i at %s with sample rate %.1f Hz...' % (obsID, tStart, sample_rate))
    idf.reset()
    
    t0 = idf.get_info('start_time')
    tDiff = tStart - datetime.utcfromtimestamp(t0)
    offset = idf.offset( tDiff.total_seconds() )
    t0 = idf.get_info('start_time')
    srate = idf.get_info('sample_rate')
    while datetime.utcfromtimestamp(t0) < tStart or srate != sample_rate:
        offset = idf.offset( 512./sample_rate )
        t0 = idf.get_info('start_time')
        srate = idf.get_info('sample_rate')
        
    print('... Found #%i at %s with sample rate %.1f Hz' % (obsID, datetime.utcfromtimestamp(t0), srate))
    tDiff = datetime.utcfromtimestamp(t0) - tStart
    duration = duration - max([0, tDiff.total_seconds()])
    
    # Number of remaining chunks (and the correction to the number of
    # frames to read in).
    nChunks = int(round(duration / args.average))
    if nChunks == 0:
        nChunks = 1
        
    # Date & Central Frequency
    beginDate = t0.datetime
    central_freq1 = idf.get_info('freq1')
    freq = numpy.fft.fftshift(numpy.fft.fftfreq(LFFT, d=1/srate))
    
    for t in range(len(antennas)//2):
        dataSets['obs%i-freq%i' % (obsID, t+1)][:] = freq + central_freq1
        
    obs = dataSets['obs%i' % obsID]
    obs.attrs['tInt'] = args.average
    obs.attrs['tInt_Unit'] = 's'
    obs.attrs['LFFT'] = LFFT
    obs.attrs['nChan'] = LFFT
    obs.attrs['RBW'] = freq[1]-freq[0]
    obs.attrs['RBW_Units'] = 'Hz'
    
    # Create the progress bar so that we can keep up with the conversion.
    pbar = progress.ProgressBarPlus(max=nChunks)
    
    data_products = ['XX', 'YY']
    done = False
    for i in xrange(nChunks):
        # Inner loop that actually reads the frames into the data array
        tInt, cTime, data = idf.read(args.average)
        if i == 0:
            print("Actual integration time is %.1f ms" % (tInt*1000.0,))
            
        # Save out some easy stuff
        dataSets['obs%i-time' % obsID][i] = (cTime[0], cTime[1])
        
        if (not args.without_sats):
            sats = ((data.real**2 + data.imag**2) >= 127**2).sum(axis=0)
            for t in range(len(antennas)//2):
                dataSets['obs%i-Saturation%i' % (obsID, t+1)][i,:] = sats[2*t+0:2*t+2]
        else:
            for t in range(len(antennas)//2):
                dataSets['obs%i-Saturation%i' % (obsID, t+1)][i,:] = -1
                
        # Calculate the spectra for this block of data and then weight the results by 
        # the total number of frames read.  This is needed to keep the averages correct.
        freq, tempSpec1 = fxc.SpecMaster(data, LFFT=LFFT, window=args.window, verbose=args.verbose, sample_rate=srate, clip_level=clip1)
        
        l = 0
        for t in range(len(antennas)//2):
            for p in data_products:
                dataSets['obs%i-%s%i' % (obsID, p, t+1)][i,:] = tempSpec1[l,:]
                l += 1
                
        # We don't really need the data array anymore, so delete it
        del(data)
        
        # Are we done yet?
        if done:
            break
            
        ## Update the progress bar and remaining time estimate
        pbar.inc()
        sys.stdout.write('%s\r' % pbar.show())
        sys.stdout.flush()
        
    pbar.amount = pbar.max
    sys.stdout.write('%s\n' % pbar.show())
    sys.stdout.flush()
    
    return True


def process_data_to_stokes(idf, antennas, tStart, duration, sample_rate, args, dataSets, obsID=1, clip1=0):
    """
    Process a chunk of data in a raw DRX file into Stokes parameters and 
    add the contents to an HDF5 file.
    """
    
    # Length of the FFT
    LFFT = args.fft_length
    
    # Find the start of the observation
    t0 = idf.get_info('start_time')
    
    print('Looking for #%i at %s with sample rate %.1f Hz...' % (obsID, tStart, sample_rate))
    idf.reset()
    
    t0 = idf.get_info('start_time')
    tDiff = tStart - datetime.utcfromtimestamp(t0)
    offset = idf.offset( tDiff.total_seconds() )
    t0 = idf.get_info('start_time')
    srate = idf.get_info('sample_rate')
    while datetime.utcfromtimestamp(t0) < tStart or srate != sample_rate:
        offset = idf.offset( 512./sample_rate )
        t0 = idf.get_info('start_time')
        srate = idf.get_info('sample_rate')
        
    print('... Found #%i at %s with sample rate %.1f Hz' % (obsID, datetime.utcfromtimestamp(t0), srate))
    tDiff = datetime.utcfromtimestamp(t0) - tStart
    duration = duration - max([0, tDiff.total_seconds()])
        
    # Number of remaining chunks (and the correction to the number of
    # frames to read in).
    nChunks = int(round(duration / args.average))
    if nChunks == 0:
        nChunks = 1
        
    # Number of remaining chunks (and the correction to the number of
    # frames to read in).
    nChunks = int(round(duration / args.average))
    if nChunks == 0:
        nChunks = 1
        
    # Date & Central Frequency
    beginDate = t0.datetime
    central_freq1 = idf.get_info('freq1')
    freq = numpy.fft.fftshift(numpy.fft.fftfreq(LFFT, d=1/srate))
    
    for t in range(len(antennas)//2):
        dataSets['obs%i-freq%i' % (obsID, t+1)][:] = freq + central_freq1
        
    obs = dataSets['obs%i' % obsID]
    obs.attrs['tInt'] = args.average
    obs.attrs['tInt_Unit'] = 's'
    obs.attrs['LFFT'] = LFFT
    obs.attrs['nChan'] = LFFT
    obs.attrs['RBW'] = freq[1]-freq[0]
    obs.attrs['RBW_Units'] = 'Hz'
    
    # Create the progress bar so that we can keep up with the conversion.
    pbar = progress.ProgressBarPlus(max=nChunks)
    
    data_products = ['I', 'Q', 'U', 'V']
    done = False
    for i in xrange(nChunks):
        # Inner loop that actually reads the frames into the data array
        tInt, cTime, data = idf.read(args.average)
        if i == 0:
            print("Actual integration time is %.1f ms" % (tInt*1000.0,))
            
        # Save out some easy stuff
        dataSets['obs%i-time' % obsID][i] = (cTime[0], cTime[1])
        
        if (not args.without_sats):
            sats = ((data.real**2 + data.imag**2) >= 127**2).sum(axis=0)
            for t in range(len(antennas)//2):
                dataSets['obs%i-Saturation%i' % (obsID, t+1)][i,:] = sats[2*t+0:2*t+2]
        else:
            for t in range(len(antennas)//2):
                dataSets['obs%i-Saturation%i' % (obsID, t+1)][i,:] = -1
                
        # Calculate the spectra for this block of data and then weight the results by 
        # the total number of frames read.  This is needed to keep the averages correct.
        freq, tempSpec1 = fxc.StokesMaster(data, antennas, LFFT=LFFT, window=args.window, verbose=args.verbose, sample_rate=srate, clip_level=clip1)
        
        for t in range(len(antennas)//2):
            for l,p in enumerate(data_products):
                dataSets['obs%i-%s%i' % (obsID, p, t+1)][i,:] = tempSpec1[l,t-1,:]
                
        # We don't really need the data array anymore, so delete it
        del(data)
        
        # Are we done yet?
        if done:
            break
            
        ## Update the progress bar and remaining time estimate
        pbar.inc()
        sys.stdout.write('%s\r' % pbar.show())
        sys.stdout.flush()
        
    pbar.amount = pbar.max
    sys.stdout.write('%s\n' % pbar.show())
    sys.stdout.flush()
    
    return True


def main(args):
    # Length of the FFT and the window to use
    LFFT = args.fft_length
    if args.bartlett:
        window = numpy.bartlett
    elif args.blackman:
        window = numpy.blackman
    elif args.hanning:
        window = numpy.hanning
    else:
        window = fxc.null_window
    args.window = window
    
    # Open the file and find good data
    idf = LWA1DataFile(args.filename, ignore_timetag_errors=args.ignore_time_errors)
    
    # Metadata
    nFramesFile = idf.get_info('nframe')
    srate = idf.get_info('sample_rate')
    antpols = idf.get_info('nantenna')
    
    # Number of frames to integrate over
    nFramesAvg = int(args.average * srate / 512) * antpols
    nFramesAvg = int(1.0 * (nFramesAvg // antpols)*512/float(LFFT))*LFFT/512 * antpols
    args.average = 1.0 * (nFramesAvg // antpols) * 512 / srate
    maxFrames = nFramesAvg
    
    # Offset into the file, if needed
    offset = idf.offset(args.skip)
    
    # Number of remaining chunks (and the correction to the number of
    # frames to read in).
    if args.metadata is not None:
        args.duration = 0
    if args.duration == 0:
        args.duration = 1.0 * nFramesFile / antpols * 512 / srate
        args.duration -= args.skip
    else:
        args.duration = int(round(args.duration * srate * antpols / 512) // antpols * 512 // srate)
    nChunks = int(round(args.duration / args.average))
    if nChunks == 0:
        nChunks = 1
    nFrames = nFramesAvg*nChunks
    
    # Date & Central Frequency
    t1  = idf.get_info('start_time')
    beginDate = t1.datetime
    central_freq1 = idf.get_info('freq1')
    
    # File summary
    print("Filename: %s" % args.filename)
    print("Date of First Frame: %s" % str(beginDate))
    print("Antenna/Pols: %i" % antpols)
    print("Sample Rate: %i Hz" % srate)
    print("Tuning Frequency: %.3f Hz" % (central_freq1,))
    print("Frames: %i (%.3f s)" % (nFramesFile, 1.0 * nFramesFile / antpols * 512 / srate))
    print("---")
    print("Offset: %.3f s (%i frames)" % (args.skip, offset))
    print("Integration: %.3f s (%i frames; %i frames per antenna/pol)" % (args.average, nFramesAvg, nFramesAvg // antpols))
    print("Duration: %.3f s (%i frames; %i frames per antenna/pol)" % (args.average*nChunks, nFrames, nFrames // antpols))
    print("Chunks: %i" % nChunks)
    print(" ")
    
    # Estimate clip level (if needed)
    if args.estimate_clip_level:
        estimate = idf.estimate_levels(sigma=5.0)
        clip1 = 1.0*sum(estimate) / len(estimate)
    else:
        clip1 = args.clip_level
        
    # Get the antennas for Stokes calculation
    if args.metadata is not None:
        try:
            project = metabundle.get_sdf(args.metadata)
            station = stations.lwa1
        except Exception as e:
            project = metabundleADP.get_sdf(args.metadata)
            station = stations.lwasv
    elif args.lwasv:
        station = stations.lwasv
    else:
        station = stations.lwa1
    antennas = station.antennas
    
    # Setup the output file
    outname = os.path.split(args.filename)[1]
    outname = os.path.splitext(outname)[0]
    outname = '%s-tbn-waterfall.hdf5' % outname
    
    if os.path.exists(outname):
        if not args.force:
            yn = raw_input("WARNING: '%s' exists, overwrite? [Y/n] " % outname)
        else:
            yn = 'y'
            
        if yn not in ('n', 'N'):
            os.unlink(outname)
        else:
            raise RuntimeError("Output file '%s' already exists" % outname)
            
    f = hdfData.create_new_file(outname)
    
    # Look at the metadata and come up with a list of observations.  If 
    # there are no metadata, create a single "observation" that covers the
    # whole file.
    obsList = {}
    if args.metadata is not None:
        try:
            project = metabundle.get_sdf(args.metadata)
        except Exception as e:
            project = metabundleADP.get_sdf(args.metadata)
            
        sdfBeam  = project.sessions[0].drx_beam
        if sdfBeam != 5:
            raise RuntimeError("Metadata is for beam #%i, but data is from beam #%i" % (sdfBeam, 5))
            
        for i,obs in enumerate(project.sessions[0].observations):
            sdfStart = mcs.mjdmpm_to_datetime(obs.mjd, obs.mpm)
            sdfStop  = mcs.mjdmpm_to_datetime(obs.mjd, obs.mpm + obs.dur)
            obsDur   = obs.dur/1000.0
            obsSR    = tbn.FILTER_CODES[obs.filter]
            
            obsList[i+1] = (sdfStart, sdfStop, obsDur, obsSR)
            
        print("Observations:")
        for i in sorted(obsList.keys()):
            obs = obsList[i]
            print(" #%i: %s to %s (%.3f s) at %.3f MHz" % (i, obs[0], obs[1], obs[2], obs[3]/1e6))
        print(" ")
            
        hdfData.fill_from_metabundle(f, args.metadata)
        
    elif args.sdf is not None:
        try:
            project = sdf.parse_sdf(args.sdf)
        except Exception as e:
            project = sdfADP.parse_sdf(args.sdf)
            
        sdfBeam  = project.sessions[0].drx_beam
        if sdfBeam != 5:
            raise RuntimeError("Metadata is for beam #%i, but data is from beam #%i" % (sdfBeam, 5))
            
        for i,obs in enumerate(project.sessions[0].observations):
            sdfStart = mcs.mjdmpm_to_datetime(obs.mjd, obs.mpm)
            sdfStop  = mcs.mjdmpm_to_datetime(obs.mjd, obs.mpm + obs.dur)
            obsDur   = obs.dur/1000.0
            obsSR    = tbn.FILTER_CODES[obs.filter]
            
            obsList[i+1] = (sdfStart, sdfStop, obsDur, obsSR)
            
        site = 'lwa1'
        if args.lwasv:
            site = 'lwasv'
        hdfData.fill_from_sdf(f, args.sdf, station=site)
        
    else:
        obsList[1] = (datetime.utcfromtimestamp(t1), datetime(2222,12,31,23,59,59), args.duration, srate)
        
        site = 'lwa1'
        if args.lwasv:
            site = 'lwasv'
        hdfData.fill_minimum(f, 1, 5, srate, station=site)
        
    if (not args.stokes):
        data_products = ['XX', 'YY']
    else:
        data_products = ['I', 'Q', 'U', 'V']
        
    for o in sorted(obsList.keys()):
        for t in range(len(antennas)//2):
            hdfData.create_observation_set(f, o, t+1, numpy.arange(LFFT, dtype=numpy.float64), int(round(obsList[o][2]/args.average)), data_products)
            
    f.attrs['FileGenerator'] = 'tbnWaterfall.py'
    f.attrs['InputData'] = os.path.basename(args.filename)
    
    # Create the various HDF group holders
    ds = {}
    for o in sorted(obsList.keys()):
        obs = hdfData.get_observation_set(f, o)
        
        ds['obs%i' % o] = obs
        ds['obs%i-time' % o] = hdfData.get_time(f, o)
        
        for t in range(len(antennas)//2):
            ds['obs%i-freq%i' % (o, t+1)] = hdfData.get_data_set(f, o, t+1, 'freq')
            for p in data_products:
                ds["obs%i-%s%i" % (o, p, t+1)] = hdfData.get_data_set(f, o, t+1, p)
            ds['obs%i-Saturation%i' % (o, t+1)] = hdfData.get_data_set(f, o, t+1, 'Saturation')
            
    # Load in the correct analysis function
    if (not args.stokes):
        process_data = process_data_to_linear
    else:
        process_data = process_data_to_stokes
        
    # Go!
    for o in sorted(obsList.keys()):
        try:
            process_data(idf, antennas, obsList[o][0], obsList[o][2], obsList[o][3], args, ds, obsID=o, clip1=clip1)
        except RuntimeError as e:
            print("Observation #%i: %s, abandoning this observation" % (o, str(e)))

    # Save the output to a HDF5 file
    f.close()
    
    # Close out the data file
    idf.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
            description='read in a TBN file and create a collection of time-averaged spectra stored as an HDF5 file called <filename>-tbn-waterfall.hdf5', 
            epilog='NOTE:  Both the -m/--metadata and -i/--sdf options provide the same additional observation information to %(prog)s so only one needs to be provided.  Specifying the -m/--metadata or -i/--sdf optiosn overrides the -d/--duration setting and the entire file is reduced.', 
            formatter_class=argparse.ArgumentDefaultsHelpFormatter
            )
    parser.add_argument('filename', type=str, 
                        help='filename to process')
    wgroup = parser.add_mutually_exclusive_group(required=False)
    wgroup.add_argument('-t', '--bartlett', action='store_true', 
                        help='apply a Bartlett window to the data')
    wgroup.add_argument('-b', '--blackman', action='store_true', 
                        help='apply a Blackman window to the data')
    wgroup.add_argument('-n', '--hanning', action='store_true', 
                        help='apply a Hanning window to the data')
    parser.add_argument('-s', '--skip', type=aph.positive_or_zero_float, default=0.0, 
                        help='skip the specified number of seconds at the beginning of the file')
    parser.add_argument('-a', '--average', type=aph.positive_float, default=1.0, 
                        help='number of seconds of data to average for spectra')
    parser.add_argument('-d', '--duration', type=aph.positive_or_zero_float, default=0.0, 
                        help='number of seconds to calculate the waterfall for; 0 for everything') 
    parser.add_argument('-q', '--quiet', dest='verbose', action='store_false',
                        help='run %(prog)s in silent mode')
    parser.add_argument('-l', '--fft-length', type=aph.positive_int, default=4096, 
                        help='set FFT length')
    parser.add_argument('-c', '--clip-level', type=aph.positive_or_zero_int, default=0,  
                        help='FFT blanking clipping level in counts; 0 disables')
    parser.add_argument('-e', '--estimate-clip-level', action='store_true', 
                        help='use robust statistics to estimate an appropriate clip level; overrides -c/--clip-level')
    parser.add_argument('-m', '--metadata', type=str, 
                        help='metadata tarball for additional information')
    parser.add_argument('-i', '--sdf', type=str, 
                        help='SDF for additional information')
    sgroup = parser.add_mutually_exclusive_group(required=False)
    sgroup.add_argument('-1', '--lwa1', action='store_true', default=True, 
                        help='data is from LWA-1; needed for -i/--sdf or when no metadata is provided')
    sgroup.add_argument('-v', '--lwasv', action='store_true', 
                        help='data is from LWA-SV; needed for -i/--sdf or when no metadata is provided')
    parser.add_argument('-f', '--force', action='store_true', 
                        help='force overwritting of existing HDF5 files')
    parser.add_argument('-k', '--stokes', action='store_true', 
                        help='generate Stokes parameters instead of XX and YY')
    parser.add_argument('-w', '--without-sats', action='store_true',
                        help='do not generate saturation counts')
    parser.add_argument('-g', '--ignore-time-errors', action='store_true', 
                        help='ignore timetag errors in the file')
    args = parser.parse_args()
    main(args)
    
