#!/usr/bin/env python

from __future__ import print_function

"""
Basic script to take a COR file, apply the station delay model, and write the
data out as a CASA measurement set.

$Rev$
$Last$
$LastChangedDate$
"""

import os
import sys
import aipy
import time
import ephem
import numpy
import getopt
from datetime import datetime

from lsl.common import stations
from lsl.reader.ldp import CORFile
from lsl.common.constants import c as speedOfLight
from lsl import astro
from lsl.imaging import utils
from lsl.correlator import uvUtils
from lsl.sim import vis as simVis
from lsl.writer import measurementset

from matplotlib import pyplot as plt


def usage(exitCode=None):
    print("""cor2MS.py - Given a COR files created by ADP, convert the file into a CASA
meaurement set.

Usage: corToMS.py [OPTIONS] file

Options:
-h, --help                  Display this help information
-s, --skip                  Skip the specified number of seconds at the beginning
                            of the file (default = 0)
-d, --duration              Number of seconds to write out data for (default = 0; 
                            run the entire file)
-o, --output                Write the combined file to the provided filename
                            (Default = auto-deterine the filename)
""")
    
    if exitCode is not None:
        sys.exit(exitCode)
    else:
        return True


def parseOptions(args):
    config = {}
    # Command line flags - default values
    config['offset'] = 0.0
    config['duration'] = 0.0
    config['output'] = None
    config['args'] = []
    
    # Read in and process the command line flags
    try:
        opts, args = getopt.getopt(args, "hs:d:o:", ["help", "skip=", "duration=", "output="])
    except getopt.GetoptError, err:
        # Print help information and exit:
        print(str(err)) # will print something like "option -a not recognized"
        usage(exitCode=2)
        
    # Work through opts
    for opt, value in opts:
        if opt in ('-h', '--help'):
            usage(exitCode=0)
        elif opt in ('-s', '--skip'):
            config['offset'] = float(value)
        elif opt in ('-d', '--duration'):
            config['duration'] = float(value)
        elif opt in ('-o', '--output'):
            config['output'] = value
        else:
            assert False
            
    # Add in arguments
    config['args'] = args
    
    # Validate
    if len(config['args']) != 1:
        raise RuntimeError("Must provide at a single file to convert")
        
    
    # Return configuration
    return config


def main(args):
    # Parse the command line
    config = parseOptions(args)
    filename = config['args'][0]
    
    # Setup the site information
    station = stations.lwasv
    ants = station.getAntennas()
    nAnt = len([a for a in ants if a.pol == 0])
    
    # Open the file and get ready to go
    idf = CORFile(filename)
    nBL = idf.getInfo('nBaseline')
    nChan = idf.getInfo('nChan')
    tInt = idf.getInfo('tInt')
    nFpO = nBL * nChan / 72
    nInts = idf.getInfo('nFrames') / nFpO
    
    jd = astro.unix_to_utcjd(idf.getInfo('tStart'))
    date = str(ephem.Date(jd - astro.DJD_OFFSET))
    centralFreq = idf.getInfo('freq1')
    centralFreq = centralFreq[len(centralFreq)/2]
    
    print("Data type:  %s" % type(idf))
    print("Samples per observations: %i" % (nFpO,))
    print("Integration Time: %.3f s" % tInt)
    print("Tuning frequency: %.3f Hz" % centralFreq)
    print("Captures in file: %i (%.1f s)" % (nInts, nInts*tInt))
    print("==")
    print("Station: %s" % station.name)
    print("Date observed: %s" % date)
    print("Julian day: %.5f" % jd)
    print(" ")
    
    # Offset into the file
    offset = idf.offset(config['offset'])
    if offset != 0.0:
        print("Skipped %.3f s into the file" % offset)
        
    # Open the file and go
    nFiles = int(config['duration'] /  tInt)
    if nFiles == 0:
        nFiles = numpy.inf
    
    fileCount = 0
    while fileCount < nFiles:
        try:
            tInt, tStart, data = idf.read(tInt)
        except Exception as e:
            print("ERROR: %s" % str(e))
            break
            
        freqs = idf.getInfo('freq1')
        beginJD = astro.unix_to_utcjd( tStart )
        beginTime = datetime.utcfromtimestamp( tStart )
        
        try:
            phase
        except NameError:
            print("Updating phasing for %.3f to %.3f MHz" % (freqs[0]/1e6, freqs[-1]/1e6))
            
            k = 0
            phase = numpy.zeros((nBL, nChan, 2, 2), dtype=numpy.complex64)
            gaix = [a.cable.gain(freqs) for a in ants if a.pol == 0]
            gaiy = [a.cable.gain(freqs) for a in ants if a.pol == 1]
            dlyx = [a.cable.delay(freqs) - a.stand.z / speedOfLight for a in ants if a.pol == 0]
            dlyy = [a.cable.delay(freqs) - a.stand.z / speedOfLight for a in ants if a.pol == 1]
            for i in xrange(nAnt):
                for j in xrange(i, nAnt):
                    phase[k,:,0,0] = numpy.exp(2j*numpy.pi*freqs*(dlyx[i] - dlyx[j])) \
                                        / numpy.sqrt(gaix[i]*gaix[j])
                    phase[k,:,0,1] = numpy.exp(2j*numpy.pi*freqs*(dlyx[i] - dlyy[j])) \
                                        / numpy.sqrt(gaix[i]*gaiy[j])
                    phase[k,:,1,0] = numpy.exp(2j*numpy.pi*freqs*(dlyy[i] - dlyx[j])) \
                                        / numpy.sqrt(gaiy[i]*gaix[j])
                    phase[k,:,1,1] = numpy.exp(2j*numpy.pi*freqs*(dlyy[i] - dlyy[j])) \
                                        / numpy.sqrt(gaiy[i]*gaiy[j])
                    
                    k += 1
                    
        for i in xrange(data.shape[-1]):
            data[...,i] *= phase
            
        # Convert to a dataDict
        try:
            blList
        except NameError:
            blList = uvUtils.getBaselines(ants[0::2], IncludeAuto=True)
            
        if config['output'] is None:
            outname = os.path.basename(filename)
            outname = os.path.splitext(outname)[0]
            outname = "%s_%i_%s.ms" % (outname, int(beginJD-astro.MJD_OFFSET), beginTime.strftime("%H_%M_%S"))
        else:
            base, ext = os.path.splitext(config['output'])
            if ext == '':
                ext = '.ms'
            outname = "%s_%i_%s%s" % (base, int(beginJD-astro.MJD_OFFSET), beginTime.strftime("%H_%M_%S"), ext)
            
        fits = measurementset.MS(outname, refTime=tStart)
        fits.setStokes(['xx', 'xy', 'yx', 'yy'])
        fits.setFrequency(freqs)
        fits.setGeometry(station, ants[0::2])
        
        obsTime = astro.unix_to_taimjd(tStart)
        fits.addDataSet(obsTime, tInt, blList, data[:,:,0,0,0], pol='xx')
        fits.addDataSet(obsTime, tInt, blList, data[:,:,0,1,0], pol='xy')
        fits.addDataSet(obsTime, tInt, blList, data[:,:,1,0,0], pol='yx')
        fits.addDataSet(obsTime, tInt, blList, data[:,:,1,1,0], pol='yy')
        fits.write()
        fits.close()
        fileCount += 1
        
    idf.close()


if __name__ == "__main__":
    main(sys.argv[1:])
    
