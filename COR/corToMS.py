#!/usr/bin/env python

"""
Basic script to take a COR file, apply the station delay model, and write the
data out as a CASA measurement set.
"""

# Python3 compatiability
from __future__ import print_function, division
import sys
if sys.version_info > (3,):
    xrange = range
    
import os
import sys
import aipy
import time
import numpy
import argparse
from datetime import datetime
from astropy.constants import c as speedOfLight
speedOfLight = speedOfLight.to('m/s').value

from lsl.common import stations
from lsl.reader.ldp import CORFile
from lsl import astro
from lsl.imaging import utils
from lsl.correlator import uvutils
from lsl.sim import vis as simVis
from lsl.writer import measurementset
from lsl.misc import parser as aph

from matplotlib import pyplot as plt




def main(args):
    # Setup the site information
    station = stations.lwasv
    ants = station.antennas
    nAnt = len([a for a in ants if a.pol == 0])
    
    phase_freq_range = [0, 0]
    for filename in args.filename:
        # Open the file and get ready to go
        idf = CORFile(filename)
        nBL = idf.get_info('nbaseline')
        nchan = idf.get_info('nchan')
        tInt = idf.get_info('tint')
        nFpO = nBL * nchan / 72
        nInts = idf.get_info('nframe') / nFpO
        
        jd = astro.unix_to_utcjd(idf.get_info('start_time'))
        date = idf.get_info('start_time').datetime
        central_freq = idf.get_info('freq1')
        central_freq = central_freq[len(central_freq)/2]
        
        print("Data type:  %s" % type(idf))
        print("Samples per observations: %i" % (nFpO,))
        print("Integration Time: %.3f s" % tInt)
        print("Tuning frequency: %.3f Hz" % central_freq)
        print("Captures in file: %i (%.1f s)" % (nInts, nInts*tInt))
        print("==")
        print("Station: %s" % station.name)
        print("Date observed: %s" % date)
        print("Julian day: %.5f" % jd)
        print(" ")
        
        # Offset into the file
        offset = idf.offset(args.skip)
        if offset != 0.0:
            print("Skipped %.3f s into the file" % offset)
            
        # Open the file and go
        nFiles = int(args.duration /  tInt)
        if nFiles == 0:
            nFiles = numpy.inf
            
        fileCount = 0
        while fileCount < nFiles:
            try:
                tInt, tStart, data = idf.read(tInt)
            except Exception as e:
                print("ERROR: %s" % str(e))
                break
                
            freqs = idf.get_info('freq1')
            beginJD = astro.unix_to_utcjd( tStart )
            beginTime = datetime.utcfromtimestamp( tStart )
            
            if freqs[0] != phase_freq_range[0] or freqs[-1] != phase_freq_range[1]:
                print("Updating phasing for %.3f to %.3f MHz" % (freqs[0]/1e6, freqs[-1]/1e6))
                phase_freq_range[0] = freqs[ 0]
                phase_freq_range[1] = freqs[-1]
                
                k = 0
                phase = numpy.zeros((nBL, nchan, 2, 2), dtype=numpy.complex64)
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
                blList  # pylint: disable=used-before-assignment
            except NameError:
                blList = uvutils.get_baselines(ants[0::2], include_auto=True)
                
            if args.output is None:
                outname = os.path.basename(filename)
                outname = os.path.splitext(outname)[0]
                outname = "%s_%i_%s.ms" % (outname, int(beginJD-astro.MJD_OFFSET), beginTime.strftime("%H_%M_%S"))
            else:
                base, ext = os.path.splitext(args.output)
                if ext == '':
                    ext = '.ms'
                outname = "%s_%i_%s%s" % (base, int(beginJD-astro.MJD_OFFSET), beginTime.strftime("%H_%M_%S"), ext)
                
            fits = measurementset.Ms(outname, ref_time=tStart)
            fits.set_stokes(['xx', 'xy', 'yx', 'yy'])
            fits.set_frequency(freqs)
            fits.set_geometry(station, ants[0::2])
            
            obsTime = astro.unix_to_taimjd(tStart)
            fits.add_data_set(obsTime, tInt, blList, data[:,:,0,0,0], pol='xx')
            fits.add_data_set(obsTime, tInt, blList, data[:,:,0,1,0], pol='xy')
            fits.add_data_set(obsTime, tInt, blList, data[:,:,1,0,0], pol='yx')
            fits.add_data_set(obsTime, tInt, blList, data[:,:,1,1,0], pol='yy')
            fits.write()
            fits.close()
            fileCount += 1
            
        idf.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='given a COR files created by ADP, convert the file into a CASA meaurement set', 
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
    parser.add_argument('filename', type=str, nargs='+', 
                        help='filename to convert')
    parser.add_argument('-s', '--skip', type=aph.positive_or_zero_float, default=0.0, 
                        help='skip the specified number of seconds at the beginning of the file')
    parser.add_argument('-d', '--duration', type=aph.positive_or_zero_float, default=0.0, 
                        help='number of seconds to write out data for')
    parser.add_argument('-o', '--output', type=str, 
                        help='write the combined file to the provided filename, auto-determine if not provided')
    args = parser.parse_args()
    main(args)
    
