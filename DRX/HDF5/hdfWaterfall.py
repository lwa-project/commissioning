#!/usr/bin/env python

"""
Given a DRX file, plot the time averaged spectra for each beam output over some 
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
import ephem
import argparse
from datetime import datetime

from lsl.reader import drx, drspec, errors
from lsl.reader.ldp import LWA1DataFile
import lsl.correlator.fx as fxc
from lsl.astro import unix_to_utcjd, DJD_OFFSET
from lsl.common import progress, stations
from lsl.common import mcs, sdf, metabundle
try:
    from lsl.common import sdfADP, metabundleADP
    adpReady = True
except ImportError:
    adpReady = False
from lsl.misc import parser as aph

import matplotlib.pyplot as plt

import data as hdfData


def bestFreqUnits(freq):
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


def processDataBatchLinear(idf, antennas, tStart, duration, sample_rate, args, dataSets, obsID=1, clip1=0, clip2=0):
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
        offset = idf.offset( 4096./sample_rate )
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
    beginDate = ephem.Date(unix_to_utcjd(t0) - DJD_OFFSET)
    central_freq1 = idf.get_info('freq1')
    central_freq2 = idf.get_info('freq2')
    freq = numpy.fft.fftshift(numpy.fft.fftfreq(LFFT, d=1/srate))
    
    dataSets['obs%i-freq1' % obsID][:] = freq + central_freq1
    dataSets['obs%i-freq2' % obsID][:] = freq + central_freq2
    
    obs = dataSets['obs%i' % obsID]
    obs.attrs['tInt'] = args.average
    obs.attrs['tInt_Unit'] = 's'
    obs.attrs['LFFT'] = LFFT
    obs.attrs['nChan'] = LFFT
    obs.attrs['RBW'] = freq[1]-freq[0]
    obs.attrs['RBW_Units'] = 'Hz'
    
    data_products = ['XX', 'YY']
    done = False
    for i in xrange(nChunks):
        # Inner loop that actually reads the frames into the data array
        print("Working on chunk %i, %i chunks remaning" % (i+1, nChunks-i-1))
        print("Working on %.1f ms of data" % (args.average*1000.0,))
        
        tInt, cTime, data = idf.read(args.average)
        if i == 0:
            print("Actual integration time is %.1f ms" % (tInt*1000.0,))
            
        # Save out some easy stuff
        dataSets['obs%i-time' % obsID][i] = float(cTime)
        
        if (not args.without_sats):
            sats = ((data.real**2 + data.imag**2) >= 49).sum(axis=1)
            dataSets['obs%i-Saturation1' % obsID][i,:] = sats[0:2]
            dataSets['obs%i-Saturation2' % obsID][i,:] = sats[2:4]
        else:
            dataSets['obs%i-Saturation1' % obsID][i,:] = -1
            dataSets['obs%i-Saturation2' % obsID][i,:] = -1
            
        # Calculate the spectra for this block of data and then weight the results by 
        # the total number of frames read.  This is needed to keep the averages correct.
        if clip1 == clip2:
            freq, tempSpec1 = fxc.SpecMaster(data, LFFT=LFFT, window=args.window, verbose=args.verbose, sample_rate=srate, clip_level=clip1)
            
            l = 0
            for t in (1,2):
                for p in data_products:
                    dataSets['obs%i-%s%i' % (obsID, p, t)][i,:] = tempSpec1[l,:]
                    l += 1
                    
        else:
            freq, tempSpec1 = fxc.SpecMaster(data[:2,:], LFFT=LFFT, window=args.window, verbose=args.verbose, sample_rate=srate, clip_level=clip1)
            freq, tempSpec2 = fxc.SpecMaster(data[2:,:], LFFT=LFFT, window=args.window, verbose=args.verbose, sample_rate=srate, clip_level=clip2)
            
            for l,p in enumerate(data_products):
                dataSets['obs%i-%s%i' % (obsID, p, 1)][i,:] = tempSpec1[l,:]
                dataSets['obs%i-%s%i' % (obsID, p, 2)][i,:] = tempSpec2[l,:]
                
        # We don't really need the data array anymore, so delete it
        del(data)
        
        # Are we done yet?
        if done:
            break
            
    return True


def processDataBatchStokes(idf, antennas, tStart, duration, sample_rate, args, dataSets, obsID=1, clip1=0, clip2=0):
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
        offset = idf.offset( 4096./sample_rate )
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
    beginDate = ephem.Date(unix_to_utcjd(t0) - DJD_OFFSET)
    central_freq1 = idf.get_info('freq1')
    central_freq2 = idf.get_info('freq2')
    freq = numpy.fft.fftshift(numpy.fft.fftfreq(LFFT, d=1/srate))
    
    dataSets['obs%i-freq1' % obsID][:] = freq + central_freq1
    dataSets['obs%i-freq2' % obsID][:] = freq + central_freq2
    
    obs = dataSets['obs%i' % obsID]
    obs.attrs['tInt'] = args.average
    obs.attrs['tInt_Unit'] = 's'
    obs.attrs['LFFT'] = LFFT
    obs.attrs['nChan'] = LFFT
    obs.attrs['RBW'] = freq[1]-freq[0]
    obs.attrs['RBW_Units'] = 'Hz'
    
    data_products = ['I', 'Q', 'U', 'V']
    done = False
    for i in xrange(nChunks):
        # Inner loop that actually reads the frames into the data array
        print("Working on chunk %i, %i chunks remaning" % (i+1, nChunks-i-1))
        print("Working on %.1f ms of data" % (args.average*1000.0,))
        
        tInt, cTime, data = idf.read(args.average)
        if i == 0:
            print("Actual integration time is %.1f ms" % (tInt*1000.0,))
            
        # Save out some easy stuff
        dataSets['obs%i-time' % obsID][i] = float(cTime)
        
        if (not args.without_sats):
            sats = ((data.real**2 + data.imag**2) >= 49).sum(axis=1)
            dataSets['obs%i-Saturation1' % obsID][i,:] = sats[0:2]
            dataSets['obs%i-Saturation2' % obsID][i,:] = sats[2:4]
        else:
            dataSets['obs%i-Saturation1' % obsID][i,:] = -1
            dataSets['obs%i-Saturation2' % obsID][i,:] = -1
            
        # Calculate the spectra for this block of data and then weight the results by 
        # the total number of frames read.  This is needed to keep the averages correct.
        if clip1 == clip2:
            freq, tempSpec1 = fxc.StokesMaster(data, antennas, LFFT=LFFT, window=args.window, verbose=args.verbose, sample_rate=srate, clip_level=clip1)
            
            for t in (1,2):
                for l,p in enumerate(data_products):
                    dataSets['obs%i-%s%i' % (obsID, p, t)][i,:] = tempSpec1[l,t-1,:]
                    
        else:
            freq, tempSpec1 = fxc.StokesMaster(data[:2,:], antennas[:2], LFFT=LFFT, window=args.window, verbose=args.verbose, sample_rate=srate, clip_level=clip1)
            freq, tempSpec2 = fxc.StokesMaster(data[2:,:], antennas[2:], LFFT=LFFT, window=args.window, verbose=args.verbose, sample_rate=srate, clip_level=clip2)
            
            for l,p in enumerate(data_products):
                dataSets['obs%i-%s%i' % (obsID, p, 1)][i,:] = tempSpec1[l,0,:]
                dataSets['obs%i-%s%i' % (obsID, p, 2)][i,:] = tempSpec2[l,0,:]
                
        # We don't really need the data array anymore, so delete it
        del(data)
        
        # Are we done yet?
        if done:
            break
            
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
    
    # Open the file and find good data (not spectrometer data)
    fh = open(args.filename, "rb")
    
    try:
        for i in xrange(5):
            junkFrame = drspec.read_frame(fh)
        raise RuntimeError("ERROR: '%s' appears to be a DR spectrometer file, not a raw DRX file" % args.filename)
    except errors.SyncError:
        fh.seek(0)
        
    # Good, we seem to have a real DRX file, switch over to the LDP interface
    fh.close()
    idf = LWA1DataFile(args.filename, ignore_timetag_errors=args.ignore_time_errors)

    # Metadata
    nFramesFile = idf.get_info('nframe')
    beam = idf.get_info('beam')
    srate = idf.get_info('sample_rate')
    beampols = idf.get_info('nbeampol')
    beams = max([1, beampols // 4])
    
    # Number of frames to integrate over
    nFramesAvg = int(args.average * srate / 4096 * beampols)
    nFramesAvg = int(1.0 * nFramesAvg / beampols*4096/float(LFFT))*LFFT/4096*beampols
    args.average = 1.0 * nFramesAvg / beampols * 4096 / srate
    maxFrames = nFramesAvg
    
    # Offset into the file, if needed
    offset = idf.offset(args.skip)
    
    # Number of remaining chunks (and the correction to the number of
    # frames to read in).
    if args.metadata is not None:
        args.duration = 0
    if args.duration == 0:
        args.duration = 1.0 * nFramesFile / beampols * 4096 / srate
        args.duration -= args.skip
    else:
        args.duration = int(round(args.duration * srate * beampols / 4096) / beampols * 4096 / srate)
    nChunks = int(round(args.duration / args.average))
    if nChunks == 0:
        nChunks = 1
    nFrames = nFramesAvg*nChunks
    
    # Date & Central Frequency
    t1  = idf.get_info('start_time')
    beginDate = ephem.Date(unix_to_utcjd(t1) - DJD_OFFSET)
    central_freq1 = idf.get_info('freq1')
    central_freq2 = idf.get_info('freq2')
    
    # File summary
    print("Filename: %s" % args.filename)
    print("Date of First Frame: %s" % str(beginDate))
    print("Beams: %i" % beams)
    print("Tune/Pols: %i" % beampols)
    print("Sample Rate: %i Hz" % srate)
    print("Tuning Frequency: %.3f Hz (1); %.3f Hz (2)" % (central_freq1, central_freq2))
    print("Frames: %i (%.3f s)" % (nFramesFile, 1.0 * nFramesFile / beampols * 4096 / srate))
    print("---")
    print("Offset: %.3f s (%i frames)" % (args.skip, offset))
    print("Integration: %.3f s (%i frames; %i frames per beam/tune/pol)" % (args.average, nFramesAvg, nFramesAvg / beampols))
    print("Duration: %.3f s (%i frames; %i frames per beam/tune/pol)" % (args.average*nChunks, nFrames, nFrames / beampols))
    print("Chunks: %i" % nChunks)
    print(" ")
    
    # Estimate clip level (if needed)
    if args.estimate_clip_level:
        estimate = idf.estimate_levels(fh, sigma=5.0)
        clip1 = (estimate[0] + estimate[1]) / 2.0
        clip2 = (estimate[2] + estimate[3]) / 2.0
    else:
        clip1 = args.clip_level
        clip2 = args.clip_level
        
    # Make the pseudo-antennas for Stokes calculation
    antennas = []
    for i in xrange(4):
        if i // 2 == 0:
            newAnt = stations.Antenna(1)
        else:
            newAnt = stations.Antenna(2)
            
        if i % 2 == 0:
            newAnt.pol = 0
        else:
            newAnt.pol = 1
            
        antennas.append(newAnt)
        
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
    
    # Look at the metadata and come up with a list of observations.  If 
    # there are no metadata, create a single "observation" that covers the
    # whole file.
    obsList = {}
    if args.metadata is not None:
        try:
            project = metabundle.get_sdf(args.metadata)
        except Exception as e:
            if adpReady:
                project = metabundleADP.get_sdf(args.metadata)
            else:
                raise e
                
        sdfBeam  = project.sessions[0].drx_beam
        spcSetup = project.sessions[0].spcSetup
        if sdfBeam != beam:
            raise RuntimeError("Metadata is for beam #%i, but data is from beam #%i" % (sdfBeam, beam))
            
        for i,obs in enumerate(project.sessions[0].observations):
            sdfStart = mcs.mjdmpm_to_datetime(obs.mjd, obs.mpm)
            sdfStop  = mcs.mjdmpm_to_datetime(obs.mjd, obs.mpm + obs.dur)
            obsDur   = obs.dur/1000.0
            obsSR    = drx.FILTER_CODES[obs.filter]
            
            obsList[i+1] = (sdfStart, sdfStop, obsDur, obsSR)
            
        print("Observations:")
        for i in sorted(obsList.keys()):
            obs = obsList[i]
            print(" #%i: %s to %s (%.3f s) at %.3f MHz" % (i, obs[0], obs[1], obs[2], obs[3]/1e6))
        print(" ")
            
        hdfData.fillFromMetabundle(f, args.metadata)
        
    elif args.sdf is not None:
        try:
            project = sdf.parse_sdf(args.sdf)
        except Exception as e:
            if adpReady:
                project = sdfADP.parse_sdf(args.sdf)
            else:
                raise e
                
        sdfBeam  = project.sessions[0].drx_beam
        spcSetup = project.sessions[0].spcSetup
        if sdfBeam != beam:
            raise RuntimeError("Metadata is for beam #%i, but data is from beam #%i" % (sdfBeam, beam))
            
        for i,obs in enumerate(project.sessions[0].observations):
            sdfStart = mcs.mjdmpm_to_datetime(obs.mjd, obs.mpm)
            sdfStop  = mcs.mjdmpm_to_datetime(obs.mjd, obs.mpm + obs.dur)
            obsDur   = obs.dur/1000.0
            obsSR    = drx.FILTER_CODES[obs.filter]
            
            obsList[i+1] = (sdfStart, sdfStop, obsDur, obsSR)
            
        site = 'lwa1'
        if args.lwasv:
            site = 'lwasv'
        hdfData.fillFromSDF(f, args.sdf, station=site)
        
    else:
        obsList[1] = (datetime.utcfromtimestamp(t1), datetime(2222,12,31,23,59,59), args.duration, srate)
        
        site = 'lwa1'
        if args.lwasv:
            site = 'lwasv'
        hdfData.fillMinimum(f, 1, beam, srate, station=site)
        
    if (not args.stokes):
        data_products = ['XX', 'YY']
    else:
        data_products = ['I', 'Q', 'U', 'V']
        
    for o in sorted(obsList.keys()):
        for t in (1,2):
            hdfData.createDataSets(f, o, t, numpy.arange(LFFT, dtype=numpy.float64), int(round(obsList[o][2]/args.average)), data_products)
            
    f.attrs['FileGenerator'] = 'hdfWaterfall.py'
    f.attrs['InputData'] = os.path.basename(args.filename)
    
    # Create the various HDF group holders
    ds = {}
    for o in sorted(obsList.keys()):
        obs = hdfData.getObservationSet(f, o)
        
        ds['obs%i' % o] = obs
        ds['obs%i-time' % o] = obs.create_dataset('time', (int(round(obsList[o][2]/args.average)),), 'f8')
        
        for t in (1,2):
            ds['obs%i-freq%i' % (o, t)] = hdfData.get_data_set(f, o, t, 'freq')
            for p in data_products:
                ds["obs%i-%s%i" % (o, p, t)] = hdfData.get_data_set(f, o, t, p)
            ds['obs%i-Saturation%i' % (o, t)] = hdfData.get_data_set(f, o, t, 'Saturation')
            
    # Load in the correct analysis function
    if (not args.stokes):
        processDataBatch = processDataBatchLinear
    else:
        processDataBatch = processDataBatchStokes
        
    # Go!
    for o in sorted(obsList.keys()):
        try:
            processDataBatch(idf, antennas, obsList[o][0], obsList[o][2], obsList[o][3], args, ds, obsID=o, clip1=clip1, clip2=clip2)
        except RuntimeError as e:
            print("Observation #%i: %s, abandoning this observation" % (o, str(e)))

    # Save the output to a HDF5 file
    f.close()
    
    # Close out the data file
    idf.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
            description='read in a DRX file and create a collection of time-averaged spectra stored as an HDF5 file called <filename>-waterfall.hdf5', 
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
    
