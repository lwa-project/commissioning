#!/usr/bin/env python

"""
Given a binary TBW health check file from PASI/LASI, covnert the data into a 
.npz file that is comptaible with the output of stationMaster.py.  These files
can be used with the smGUI.py utility.
"""

# Python3 compatiability
from __future__ import print_function, division
import sys
if sys.version_info > (3,):
    xrange = range
    
import os
import re
import sys
import numpy
import struct
try:
    from urllib2 import urlopen
except ImportError:
    from urllib.request import urlopen
import argparse
import tempfile
from datetime import datetime
from xml.etree import ElementTree
from BeautifulSoup import BeautifulSoup
    

from lsl.common.stations import parse_ssmif
from lsl.common.progress import ProgressBar


def parseIndex(index):
    """
    Parse the archive listing of SSMIF version and return a list of 
    filename/date tuples.
    """
    
    # Find the table
    start = index.find('<table>')
    stop  = index.find('</table>')
    index = index[start:stop+8]
    
    # Clean it up in such a way that ElementTree can parse it
    myMassage = [(re.compile('<!([^--].*)>'), lambda match: '<!--' + match.group(1) + '-->'), 
            (re.compile('<hr>'), lambda match: ''), 
            (re.compile('&nbsp;'), lambda match: ' '), 
            (re.compile('<a.*?>(.*)</a>'), lambda mtch: mtch.group(1))]
    soup = BeautifulSoup(index, markupMassage=myMassage)
    index = soup.prettify()
    
    # Parse it
    table = ElementTree.XML(index)
    rows = iter(table)
    
    # Extract the SSMIF entries
    versions = []
    for row in rows:
        values = [col.text for col in row]
        if len(values) != 5:
            continue
            
        filename = values[1].lstrip().rstrip()
        if filename[:5] != 'SSMIF':
            continue
        if filename.find('CURRENT') != -1:
            continue
            
        date = filename.split('_', 1)[1]
        date = date.split('.', 1)[0]
        date = date.split('_', 1)[1]
        date = datetime.strptime('%s 000000' % date, '%y%m%d %H%M%S')
        versions.append( (filename, date) )
        
    # Done
    return versions


class DynamicSSMIF(object):
    def __init__(self, filename, dt):
        self.filename = filename
        self.dt = dt
        
    def __lt__(self, y):
        return self.dt < y.dt
        
    def __le__(self, y):
        return self.dt <= y.dt
        
    def __gt__(self, y):
        return self.dt > y.dt
        
    def __ge__(self, y):
        self.dt >= y.dt
        
    def __eq__(self, y):
        self.dt == y.dt
        
    def __ne__(self, y):
        self.dt != y.dt
        
    def __cmp__(self, y):
        if self.dt > y:
            return 1
        elif self.dt < y:
            return -1
        else:
            return 0
            
    def _retrieve(self):
        """
        Pull the file from the archive, parse it, and save the various results.
        """
        
        # Pull the data from the archive
        ah = urlopen("https://lda10g.alliance.unm.edu/metadata/lwasv/ssmif/%s" % self.filename)
        contents = ah.read()
        ah.close()
        
        # Save it to a file
        _, filename = tempfile.mkstemp(suffix='.txt', prefix='SSMIF')
        fh = open(filename, 'wb')
        fh.write(contents)
        fh.close()
        
        # Parse the SSMIF
        station = parse_ssmif(filename)
        
        # Cleanup
        os.unlink(filename)
        
        # Save and done
        self.contents = contents
        self.station = station
        return True
        
    def getDate(self):
        return self.dt
        
    def getContents(self):
        try:
            contents = self.contents
        except AttributeError:
            self._retrieve()
            contents = self.contents
        contents = contents.split('\n')
        
        return contents
        
    def get_station(self):
        try:
            station = self.station
        except AttributeError:
            self._retrieve()
            station = self.station
            
        return station


def loadSSMIFCache():
    """
    Populate the a list with DynamicSSMIF instances, one for each LWA1 
    SSMIF found in the archive.
    """
    
    # Retrieve the list
    ah = urlopen("https://lda10g.alliance.unm.edu/metadata/lwasv/ssmif/")
    index = ah.read()
    ah.close()
    
    # Parse
    versions = parseIndex(index)
    
    # Fill the list
    ssmifCache = []
    for i,(filename,date) in enumerate(versions):
        ssmifCache.append( DynamicSSMIF(filename, date) )
        
    # Sort and reverse
    ssmifCache.sort()
    ssmifCache.reverse()
    
    # Done
    return ssmifCache


def loadHealthCheckFile(filename):
    """
    Given a binary file created by PASI/LASI as part of a TBW health check 
    and parse the file into a three-element tuple of:
    * datetime of the capture, 
    * a numpy array of frequncies in Hz, and
        * a 2-D numpy array of the integrated spectra.
    """
    
    # Get the date/time of the capture
    date = os.path.basename(filename)
    date, time, junk = date.split('_', 2)
    dt = datetime.strptime('%s %s' % (date, time), '%y%m%d %H%M%S')
    
    # Open the file and parse out the spectra
    fh = open(filename, 'rb')
    nStands, nchans = struct.unpack('ll', fh.read(16))
    specs = numpy.fromfile(fh, count=nStands*2*nchans, dtype=numpy.float32)
    specs = specs.reshape(nStands*2, nchans)
    
    # Get the corresponding frequencies
    freqs = numpy.arange(nchans) / float(nchans - 1) * 196e6 / 2
    
    return dt, freqs, specs


def main(args):
    # Get the SSMIF cache
    ssmifCache = loadSSMIFCache()
    
    # Go!
    for filename in args.filename:
        ## Report
        print("Working on '%s'..." % os.path.basename(filename))
        
        ## Create the output filename and figure out if we need to overwrite if
        outfile = os.path.basename(filename)
        outfile = "%s.npz" % os.path.splitext(outfile)[0]
        if os.path.exists(outfile) and not args.force:
            print("  ERROR: Output file '%s' already exists, skipping" % outfile)
            continue
            
        ## Get the data from the file and create a masterSpectra array
        beginDate, freq, spec = loadHealthCheckFile(filename)
        masterSpectra = numpy.zeros((1,spec.shape[0],spec.shape[1]), dtype=spec.dtype)
        masterSpectra[0,:,:] = spec
        antpols, nchan = spec.shape
        
        ## Report on the file's contents
        if args.verbose:
            print("  Capture Date: %s" % beginDate)
            print("  Antenna/pols.: %i" % antpols)
            print("  Channels: %i" % nchan)
            
        ## Get the SSMIF that we need for this file
        found = False
        for ssmif in ssmifCache:
            if beginDate >= ssmif.getDate():
                found = True
                break
        if found:
            if args.verbose:
                print("  Using SSMIF '%s' for mappings" % ssmif.filename)
        else:
            print("  ERROR: Cannot find a suitable SSMIF for %s, skipping" % filename)
            continue
            
        ## Pull out the metadata we need
        station = ssmif.get_station()
        ssmifContents = ssmif.getContents()
        antennas = station.antennas
        
        ## Compute the 1 ms average power and the data range within each 1 ms window
        ## NOTE:  This is a dummy operation since we can't do this with the health
        ##        check data.
        subSize = 1960
        nsegments = 2*subSize / subSize
        avgPower = numpy.zeros((antpols, nsegments), dtype=numpy.float32)
        dataRange = numpy.zeros((antpols, nsegments, 3), dtype=numpy.int16)
        adcHistogram = numpy.zeros((antpols, 4096), dtype=numpy.int32)
        
        ## Apply the cable loss corrections, if requested
        if True:
            for s in xrange(masterSpectra.shape[1]):
                currGain = antennas[s].cable.gain(freq)
                for c in xrange(masterSpectra.shape[0]):
                    masterSpectra[c,s,:] /= currGain
                    
        ## Estimate the dipole resonance frequencies
        if args.verbose:
            print("  Computing dipole resonance frequencies")
        pb = ProgressBar(max=spec.shape[0])
        resFreq = numpy.zeros(spec.shape[0])
        toCompare = numpy.where( (freq>31e6) & (freq<70e6) )[0]
        for i in xrange(spec.shape[0]):
            bestOrder = 0
            bestRMS = 1e34
            for j in xrange(3, 12):
                coeff = numpy.polyfit(freq[toCompare]/1e6, numpy.log10(spec[i,toCompare])*10, j)
                fit = numpy.polyval(coeff, freq[toCompare]/1e6)
                rms = ((fit - numpy.log10(spec[i,toCompare])*10)**2).sum()
                if rms < bestRMS:
                    bestOrder = j
                    bestRMS = rms
                    
            coeff = numpy.polyfit(freq[toCompare]/1e6, numpy.log10(spec[i,toCompare])*10, bestOrder)
            fit = numpy.polyval(coeff, freq[toCompare]/1e6)
            try:
                resFreq[i] = freq[toCompare[numpy.where( fit == fit.max() )[0][0]]] / 1e6
            except:
                pass
                
            pb.inc(amount=1)
            if args.verbose and pb.amount != 0 and pb.amount % 10 == 0:
                sys.stdout.write('  '+pb.show()+'\r')
                sys.stdout.flush()
        if args.verbose:
            sys.stdout.write('  '+pb.show()+'\r')
            sys.stdout.write('\n')
            sys.stdout.flush()
            
        ## Save to disk
        numpy.savez(outfile, date=str(beginDate), freq=freq, masterSpectra=masterSpectra, resFreq=resFreq, 
                    avgPower=avgPower, dataRange=dataRange, adcHistogram=adcHistogram, ssmifContents=ssmifContents)
        if args.verbose:
            print("  Saved %.1f MB to '%s'" % (os.path.getsize(outfile)/1024.0**2, os.path.basename(outfile)))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='read in binary TBF health check file and convert it into a .npz file compatible with smGUI.py', 
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
    parser.add_argument('filename', type=str, nargs='+', 
                        help='filename to convert')
    parser.add_argument('-f', '--force', action='store_true', 
                        help='remake the NPZ file, even if it exists')
    parser.add_argument('-q', '--quiet', dest='verbose', action='store_false', 
                        help='run %(prog)s in silent mode')
    args = parser.parse_args()
    main(args)
    
