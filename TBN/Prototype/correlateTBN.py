#!/usr/bin/env python

"""
Example script that reads in TBN data from the prototype system at the north
arm and runs a cross-correlation on it.  The results are saved in the FITS 
IDI format.
"""

# Python3 compatiability
from __future__ import print_function, division
import sys
if sys.version_info > (3,):
    xrange = range
    
import os
import re
import sys
import time
import ephem
import numpy
import getopt
from datetime import datetime, timedelta, tzinfo

from lsl.reader.buffer import TBNFrameBuffer
from lsl import astro
from lsl.common import stations
from lsl.common import dp as dp_common
from lsl.statistics import robust
from lsl.reader import tbn
from lsl.reader import errors
from lsl.correlator import uvutil
from lsl.correlator import fx as fxc
from lsl.writer import fitsidi

import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter

import random
from collections import deque


class UTC(tzinfo):
    """tzinfo object for UTC time."""

    def utcoffset(self, dt):
        return timedelta(0)

    def tzname(self, dt):
        return "UTC"

    def dst(self, dt):
        return timedelta(0)


def usage(exitCode=None):
    print("""correlateTBN.py - cross-correlate data in a TBN file from the
prototype system at the north arm.

Usage: correlateTBN.py [OPTIONS] file

Options:
-h, --help             Display this help information
-m, --metadata         Name of SSMIF file to use for mappings
-f, --fft-length       Set FFT length (default = 256)
-t, --avg-time         Window to average visibilities in time (seconds; 
                    default = 6 s)
-s, --samples          Number of average visibilities to generate
                    (default = 10)
-o, --offset           Seconds to skip from the beginning of the file
-q, --quiet            Run correlateTBN in silent mode
-x, --xx               Compute only the XX polarization product (default)
-y, --yy               Compute only the YY polarization product
-2, --two-products     Compute both the XX and YY polarization products
-4, --four-products    Compute all for polariation products:  XX, YY, XY, 
                    and YX.
""")
    
    if exitCode is not None:
        sys.exit(exitCode)
    else:
        return True


def parseConfig(args):
    config = {}
    # Command line flags - default values
    config['SSMIF'] = ''
    config['avgTime'] = 6
    config['LFFT'] = 256
    config['samples'] = 10
    config['offset'] = 0
    config['verbose'] = True
    config['products'] = ['xx',]
    config['args'] = []

    # Read in and process the command line flags
    try:
        opts, arg = getopt.getopt(args, "hm:ql:t:s:o:24xy", ["help", "metadata=", "quiet", "fft-length=", "avg-time=", "samples=", "offset=", "two-products", "four-products", "xx", "yy"])
    except getopt.GetoptError as err:
        # Print help information and exit:
        print(str(err)) # will print something like "option -a not recognized"
        usage(exitCode=2)
    
    # Work through opts
    for opt, value in opts:
        if opt in ('-h', '--help'):
            usage(exitCode=0)
        elif opt in ('-m', '--metadata'):
            config['SSMIF'] = value
        elif opt in ('-q', '--quiet'):
            config['verbose'] = False
        elif opt in ('-l', '--fft-length'):
            config['LFFT'] = int(value)
        elif opt in ('-t', '--avg-time'):
            config['avgTime'] = float(value)
        elif opt in ('-s', '--samples'):
            config['samples'] = int(value)
        elif opt in ('-o', '--offset'):
            config['offset'] = int(value)
        elif opt in ('-2', '--two-products'):
            config['products'] = ['xx', 'yy']
        elif opt in ('-4', '--four-products'):
            config['products'] = ['xx', 'yy', 'xy', 'yx']
        elif opt in ('-x', '--xx'):
            config['products'] = ['xx',]
        elif opt in ('-y', '--yy'):
            config['products'] = ['yy',]
        else:
            assert False
    
    # Add in arguments
    config['args'] = arg

    # Return configuration
    return config


def processChunk(fh, site, good, filename, intTime=6.0, LFFT=64, overlap=1, central_freq=49.0e6, sample_rate=dp_common.fS, pols=['xx',], ChunkSize=300):
    """
    Given a filehandle pointing to some TBN data and various parameters for
    the cross-correlation, write cross-correlate the data and save it to a file.
    """

    # Get antennas
    antennas = []
    for a in site.antennas:
        if a.digitizer != 0:
            antennas.append(a)

    # Create the FrameBuffer instance and make sure that the frames coming out of
    # buffer.get() and buffer.flush() are in stand/polariation order
    buffer = TBNFrameBuffer(stands=range(1,20//2+1), pols=[0, 1], reorder=False)

    # Create the list of good digitizers and a digitizer to Antenna instance mapping.  
    # These are:
    #  mapper  -> mapping of digitizer number to array location
    #  mapper2 -> mapping of Antenna instance to array location
    mapper = [antennas[i].digitizer for i in good]
    mapper2 = [antennas[i] for i in good]
    
    # Find out how many frames to work with at a time.  This number is per stand/pol so
    # it is equivalent to the numper of successful buffer.get() calls that we should get.
    nFrames = int(intTime*sample_rate/512)

    # Main loop over the input file to read in the data and organize it.  Several control 
    # variables are defined for this:
    #  ref_time -> time (in seconds since the UNIX epoch) for the first data set
    #  setTime -> time (in seconds since the UNIX epoch) for the current data set
    ref_time = 0.0
    setTime = 0.0
    wallTime = time.time()
    for s in xrange(ChunkSize):
        iTime = 0
        count = [0 for i in good]
        data = numpy.zeros((len(good),nFrames*512), dtype=numpy.complex64)
        
        # Loop over append()'s/get()'s of the buffer
        i = 0
        doFlush = False
        while i < nFrames:
            # Read in the next frame and anticipate any problems that could occur
            cFrames = deque()
            for l in xrange(20):
                try:
                    cFrames.append( tbn.read_frame(fh) )
                except errors.EOFError:
                    doFlush = True
                    break
                except errors.SyncError:
                    print("WARNING: Mark 5C sync error on frame #%i" % (int(fh.tell())/tbn.FRAME_SIZE-1))
                    continue

            buffer.append(cFrames)
            cFrames = buffer.get()
            
            # Continue adding frames if nothing comes out.
            if cFrames is None:
                continue
    
            # If something comes out, add it to the data array
            for cFrame in cFrames:
                stand,pol = cFrame.header.id
                aStand = 2*(stand-1)+pol + 1
                
                if i == 0:
                    setTime = cFrame.time
                    if s == 0:
                        ref_time = setTime
                        
                try:
                    aStand = mapper.index(aStand)
                    data[aStand,  count[aStand]*512:(count[aStand]+1)*512] = cFrame.payload.data
                    count[aStand] = count[aStand] + 1
                except ValueError:
                    pass
                    
            i += 1
        
        # If we have encountered an EOF, flush the buffer
        if doFlush:
            for cFrames in buffer.flush():
                for cFrame in cFrames:
                    stand,pol = cFrame.header.id
                    aStand = 2*(stand-1)+pol + 1
                    
                    if i == 0:
                        setTime = cFrame.time
                        if s == 0:
                            ref_time = setTime
                    
                    try:
                        aStand = mapper.index(aStand)
                        # It could be possible for the buffer to contain more data than we actually 
                        # need.  Thus, we wrap the unload of the I/Q data into the array with a 
                        # try...expect block to mitigate problems.
                        try:
                            data[aStand,  count[aStand]*512:(count[aStand]+1)*512] = cFrame.payload.data
                            count[aStand] = count[aStand] + 1
                        except:
                            pass
                    except ValueError:
                        pass
                    
                i += 1

        # Setup the set time as a python datetime instance so that it can be easily printed
        setDT = datetime.utcfromtimestamp(setTime)
        setDT.replace(tzinfo=UTC())
        print("Working on set #%i (%.3f seconds after set #1 = %s)" % ((s+1), (setTime-ref_time), setDT.strftime("%Y/%m/%d %H:%M:%S.%f")))

        # Loop over polarization products
        for pol in pols:
            print("->  %s" % pol)
            blList, freq, vis = fxc.FXMaster(data, mapper2, LFFT=LFFT, overlap=overlap, include_auto=True, verbose=False, sample_rate=sample_rate, central_freq=central_freq, Pol=pol, return_baselines=True, gain_correct=False)

            # Select the right range of channels to save
            toUse = numpy.where( (freq>10.0e6) & (freq<93.0e6) )
            toUse = toUse[0]

            # If we are in the first polarazation product of the first iteration,  setup
            # the FITS IDI file.
            if s  == 0 and pol == pols[0]:
                pol1, pol2 = fxc.pol_to_pols(pol)
                
                fits = fitsidi.Idi(filename, ref_time=ref_time)
                fits.set_stokes(pols)
                fits.set_frequency(freq[toUse])
                fits.set_geometry(site, [a for a in mapper2 if a.pol == pol1])

            # Convert the setTime to a MJD and save the visibilities to the FITS IDI file
            obsTime = astro.unix_to_taimjd(setTime)
            fits.add_data_set(obsTime, 512*nFrames/sample_rate, blList, vis[:,toUse], pol=pol)
        print("->  Cummulative Wall Time: %.3f s (%.3f s per integration)" % ((time.time()-wallTime), (time.time()-wallTime)/(s+1)))

    # Cleanup after everything is done
    fits.write()
    fits.close()
    del(fits)
    del(data)
    del(vis)
    return True


def main(args):
    # Parse command line options
    config = parseConfig(args)
    filename = config['args'][0]

    # Length of the FFT
    LFFT = config['LFFT']

    # Setup the LWA station information
    if config['SSMIF'] != '':
        station = stations.parse_ssmif(config['SSMIF'])
    else:
        station = stations.lwana
    antennas = []
    for a in station.antennas:
        if a.digitizer != 0:
            antennas.append(a)
            
    fh = open(filename, "rb", buffering=tbn.FRAME_SIZE*10000)
    test = tbn.read_frame(fh)
    if not test.header.is_tbn:
        raise errors.notTBNError()
    central_freq = test.central_freq
    fh.seek(0)

    jd = astro.unix_to_utcjd(test.time.unix)
    date = str(ephem.Date(jd - astro.DJD_OFFSET))
    nFpO = len(antennas)
    sample_rate = tbn.get_sample_rate(fh)
    nInts = os.path.getsize(filename) // tbn.FRAME_SIZE // nFpO

    # Get valid stands for both polarizations
    goodX = []
    goodY = []
    for i in xrange(len(antennas)):
        ant = antennas[i]
        if ant.combined_status != 33:
            continue
        elif ant.pol == 0:
            goodX.append(ant)
        else:
            goodY.append(ant)
    
    # Now combine both lists to come up with stands that
    # are in both so we can form the cross-polarization 
    # products if we need to
    good = []
    for antX in goodX:
        for antY in goodY:
            if antX.stand.id == antY.stand.id:
                good.append( antX.digitizer-1 )
                good.append( antY.digitizer-1 )
    
    # Report on the valid stands found.  This is a little verbose,
    # but nice to see.
    print("Found %i good stands to use" % (len(good)//2,))
    for i in good:
        print("%3i, %i @ %i" % (antennas[i].stand.id, antennas[i].pol, antennas[i].digitizer))

    # Number of frames to read in at once and average
    nFrames = int(config['avgTime']*sample_rate/512)
    nSkip = int(config['offset']*sample_rate/512)
    fh.seek(nSkip*len(antennas)*tbn.FRAME_SIZE)
    nSets = os.path.getsize(filename) // tbn.FRAME_SIZE // nFpO // nFrames
    nSets = nSets - nSkip // nFrames

    print("TBN Data:  %s" % test.header.is_tbn)
    print("Samples per observations: %i per pol." % (nFpO//2))
    print("Filter code: %i" % tbn.get_sample_rate(fh, nframe=nFpO, filter_code=True))
    print("Sampling rate: %i Hz" % sample_rate)
    print("Tuning frequency: %.3f Hz" % central_freq)
    print("Captures in file: %i (%.1f s)" % (nInts, nInts*512 / sample_rate))
    print("==")
    print("Station: %s" % station.name)
    print("Date observed: %s" % date)
    print("Julian day: %.5f" % jd)
    print("Integration Time: %.3f s" % (512*nFrames/sample_rate))
    print("Number of integrations in file: %i" % nSets)

    # Make sure we don't try to do too many sets
    if config['samples'] > nSets:
        config['samples'] = nSets

    # Loop over junks of 300 integrations to make sure that we don't overflow 
    # the FITS IDI memory buffer
    s = 0
    leftToDo = config['samples']
    basename, ext = os.path.splitext(filename)
    while leftToDo > 0:
        fitsFilename = "%s.FITS_%i" % (basename, (s+1),)
        
        if leftToDo > 300:
            chunk = 300
        else:
            chunk = leftToDo
        
        processChunk(fh, station, good, fitsFilename, intTime=config['avgTime'], LFFT=config['LFFT'], 
                     overlap=1, central_freq=central_freq, sample_rate=sample_rate, 
                     pols=config['products'], ChunkSize=chunk)

        s += 1
        leftToDo = leftToDo - chunk

    fh.close()


if __name__ == "__main__":
    main(sys.argv[1:])
