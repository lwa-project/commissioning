#!/usr/bin/env python

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


def main(args):
    filename = args[0]
    
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
    
    print "Data type:  %s" % type(idf)
    print "Samples per observations: %i" % (nFpO,)
    print "Integration Time: %.3f s" % tInt
    print "Tuning frequency: %.3f Hz" % centralFreq
    print "Captures in file: %i (%.1f s)" % (nInts, nInts*tInt)
    print "=="
    print "Station: %s" % station.name
    print "Date observed: %s" % date
    print "Julian day: %.5f" % jd
    print " "
    
    print idf.offset(600), '->', str(ephem.Date(astro.unix_to_utcjd(idf.getInfo('tStart'))-astro.DJD_OFFSET))
    
    # Open the file and go
    q = 0
    while True:
        t0 = time.time()
        print 'Read'
        tInt, tStart, data = idf.read(tInt)
        q += 1
            
        t1 = time.time()
        print 'Metadata'
        freqs = idf.getInfo('freq1')
        beginJD = astro.unix_to_utcjd( tStart )
        
        t2 = time.time()
        try:
            phase
        except NameError:
            print 'Update Phasing for %.3f to %.3f MHz' % (freqs[0]/1e6, freqs[-1]/1e6)
            
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
                    
        t3 = time.time()
        print 'Phase'
        for i in xrange(data.shape[-1]):
            data[...,i] *= phase
            
        # Convert to a dataDict
        t4 = time.time()
        print 'Convert'
        try:
            blList
        except NameError:
            blList = uvUtils.getBaselines(ants[0::2], IncludeAuto=True)
            
        outname = os.path.basename(filename)
        outname = os.path.splitext(outname)[0]
        outname = "%s_%i.ms" % (outname, q)
        
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
        
        t5 = time.time()
        print t5-t0, '->', t1-t0, t2-t1, t3-t2, t4-t3, t5-t4
        
    fh.close()


if __name__ == "__main__":
    main(sys.argv[1:])
    