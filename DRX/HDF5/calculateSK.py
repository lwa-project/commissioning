#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
calculateSK.py - Read in a DRX/HDF5 wataerfall file and calculate the pseudo-
spectral kurtosis (pSK) for each observation (XX and YY or I).  The pSK values
are "pseudo" since the spectra contained within the HDF5 file are averaged 
over N FFTs where N > 1.  The pSK value for each channel is calculated by 
examining M consecutive spectra.

Note:  Although the pSK method works reasonable well for short integration times
(<~0.3 s) the method may break down for much longer integration times.

Usage:
./calculateSK.py [OPTIONS] file

$Rev$
$LastChangedBy$
$LastChangedDate$
"""

import os
import sys
import h5py
import numpy
import getopt
from datetime import datetime

from lsl.statistics import robust, kurtosis


def usage(exitCode=None):
    print """calculateSK.py - Read in DRX/HDF5 waterfall file and calculate the pseudo-
spectral kurtosis (pSK).

Usage: calculateSK.py [OPTIONS] file

Options:
-h, --help                Display this help information
-d, --duration            pSK update interval (default = 10s)
-n, --no-update           Do not add the pSK information to the HDF5 file
-g, --generate-mask       Added the pSK information to the file as a mask rather
                        than pSK values.  The masking parameters are 
                        controlled by the -t/--threshold, -f/--fill, and 
                        -m/--merge flags.
-t, --threshold           Masking threshold, in sigma, if the -g/--generate-mask 
                        flag is set (default = 4.0)
-f, --fill                Fill in the secondary polarizations (XY, YX; Q, U, V; 
                        RL, LR) using the mask from the primaries (XX, YY; I; 
                        RR, LL) (default = no)
-m, --merge               Merge the new mask table with the existing one, if it
                        exists (default = replace)
                        
Note: The -g/--generate-mask flag overrides the -n/--no-update flag.
"""
    
    if exitCode is not None:
        sys.exit(exitCode)
    else:
        return True


def parseOptions(args):
    config = {}
    # Command line flags - default values
    config['duration'] = 10.0
    config['update'] = True
    config['mask'] = False
    config['threshold'] = 4.0
    config['fill'] = False
    config['merge'] = False
    config['args'] = []
    
    # Read in and process the command line flags
    try:
        opts, args = getopt.getopt(args, "hd:ngt:fm", ["help", "duration=", "no-update", "generate-mask", "threshold=", "fill", "merge"])
    except getopt.GetoptError, err:
        # Print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        usage(exitCode=2)
        
    # Work through opts
    for opt, value in opts:
        if opt in ('-h', '--help'):
            usage(exitCode=0)
        elif opt in ('-d', '--duration'):
            config['duration'] = float(value)
        elif opt in ('-n', '--no-update'):
            config['update'] = False
        elif opt in ('-g', '--generate-mask'):
            config['update'] = True
            config['mask'] = True
        elif opt in ('-t', '--threshold'):
            config['threshold'] = float(value)
        elif opt in ('-f', '--fill'):
            config['fill'] = True
        elif opt in ('-m', '--merge'):
            config['merge'] = True
        else:
            assert False
            
    # Add in arguments
    config['args'] = args
    
    # Return configuration
    return config


def main(args):
    config = parseOptions(args)
    
    filename = config['args'][0]
    
    # Open the file and loop over the observation sections
    h = h5py.File(filename, mode='a' if config['update'] else 'r')
    for obsName in h.keys():
        obs = h.get(obsName, None)
        
        # Load in the information we need to calculate the pseudo-spectral kurtosis (pSK)
        tInt = obs.attrs['tInt']
        LFFT = obs.attrs['LFFT']
        srate = obs.attrs['sampleRate']
        
        skN = int(tInt*srate / LFFT)
        chunkSize = int(round(config['duration']/tInt))
        
        print "Staring Observation #%i" % int(obsName.replace('Observation', ''))
        print "  Sample Rate: %.1f Hz" % srate
        print "  LFFT: %i" % LFFT
        print "  tInt: %.3f s" % tInt
        print "  Integrations per spectrum: %i" % skN
        
        time = obs.get('time', None)[:]
        tuning1 = obs.get('Tuning1', None)
        tuning2 = obs.get('Tuning2', None)
        
        # Get all of the avaliable data products
        dataProducts = list(tuning1)
        for toRemove in ('Mask', 'Saturation', 'SpectralKurtosis', 'freq'):
            try:
                del dataProducts[dataProducts.index(toRemove)]
            except ValueError:
                pass
                
        # What's this?  Stokes I is the sum of XX and YY so there are 
        # effectively twice as many FFTs per spectrum relative to XX and
        # YY
        nAdjust = {'XX': 1, 'YY': 1, 'I': 2, 'RR': 1, 'LL': 1}
        
        # Loop over data products - primary
        for dp in ('XX', 'YY', 'I', 'RR', 'LL'):
            for t,tuning in enumerate((tuning1, tuning2)):
                data = tuning.get(dp, None)
                if data is None:
                    continue
                    
                # Loop over chunks for computing pSK
                sk = numpy.zeros_like(data)
                
                section = numpy.empty((chunkSize,data.shape[1]), dtype=data.dtype)
                for i in xrange(time.size/chunkSize+1):
                    start = i*chunkSize
                    stop = start + chunkSize
                    if stop >= time.size:
                        stop = time.size - 1
                    if start == stop:
                        continue
                        
                    ## Compute the pSK for each channel in both tunings
                    selection = numpy.s_[start:stop, :]
                    try:
                        data.read_direct(section, selection)
                    except TypeError:
                        section = data[start:stop,:]
                    for j in xrange(section.shape[1]):
                        sk[start:stop,j] = kurtosis.spectralPower(section[:,j], N=skN*nAdjust[dp])
                        
                # Report
                print "  => %s-%i SK Mean: %.3f" % (dp, t+1, numpy.mean(sk))
                print "     %s-%i SK Std. Dev.: %.3f" % (dp, t+1, numpy.std(sk))
                
                # Save the pSK information to the HDF5 file if we need to
                if config['update']:
                    h.attrs['FileLastUpdated'] = datetime.utcnow().strftime("UTC %Y/%m/%d %H:%M:%S")
                    
                    if config['mask']:
                        ## Calculate the expected pSK limits for the threshold
                        kLow, kHigh = kurtosis.getLimits(config['threshold'], chunkSize, N=skN*nAdjust[dp])
                        
                        ## Generate the mask arrays
                        maskTable = numpy.where( (sk < kLow) | (sk > kHigh), True, False )
                        
                        ## Report
                        print "     => %s-%i Mask Fraction: %.1f%%" % (dp, t+1, 100.0*maskTable.sum()/maskTable.size,)
                        
                        ## Pull out the Mask group from the HDF5 file
                        mask = tuning.get('Mask', None)
                        if mask is None:
                            mask = tuning.create_group('Mask')
                            for p in dataProducts:
                                mask.create_dataset(p, tuning[p].shape, 'bool')
                        maskDP = mask.get(dp, None)
                        if config['merge']:
                            maskDP[:,:] |= maskTable.astype(numpy.bool)
                        else:
                            maskDP[:,:] = maskTable.astype(numpy.bool)
                            
                    else:
                        sks = tuning.get('SpectralKurtosis', None)
                        if sks is None:
                            sks = tuning.create_group('SpectralKurtosis')
                        try:
                            sks.create_dataset(dp, sk.shape, 'f4')
                        except:
                            pass
                        sks[dp][:,:] = sk
                        sks.attrs['N'] = skN*nAdjust[dp]
                        sks.attrs['M'] = chunkSize
                        
        if config['mask'] and config['fill']:	
            # Loop over data products - secondary
            for dp in dataProducts:
                for t,tuning in enumerate((tuning1, tuning2)):
                    ## Jump over the primary polarizations
                    if dp in ('XX', 'YY', 'I', 'RR', 'LL'):
                        continue
                        
                    if dp in ('Q', 'U', 'V'):
                        ## Case 1 - Flag Stokes Q, U, and V off I
                        
                        ## Pull out the Mask group from the HDF5 file
                        mask = tuning.get('Mask', None)
                        maskDP = mask.get(dp, None)
                        maskBase = mask.get('I', None)
                        if config['merge']:
                            maskDP[:,:] |= maskBase[:,:]
                        else:
                            maskDP[:,:] = maskBase[:,:]
                            
                    elif dp in ('XY', 'YX'):
                        ## Case 2 - Flag XY and YX off XX and YY
                        
                        ## Pull out the Mask group from the HDF5 file
                        mask = tuning.get('Mask', None)
                        maskDP = mask.get(dp, None)
                        maskBase1 = mask.get('XX', None)
                        maskBase2 = mask.get('YY', None)
                        if config['merge']:
                            maskDP[:,:] |= (maskBase1[:,:] | maskBase2[:,:])
                        else:
                            maskDP[:,:] = (maskBase1[:,:] | maskBase2[:,:])
                            
                    elif dp in ('RL', 'LR'):
                        ## Case 3 - Flag RL and LR off RR and LL
                        
                        ## Pull out the Mask group from the HDF5 file
                        mask = tuning.get('Mask', None)
                        maskDP = mask.get(dp, None)
                        maskBase1 = mask.get('RR', None)
                        maskBase2 = mask.get('LL', None)
                        if config['merge']:
                            maskDP[:,:] |= (maskBase1[:,:] | maskBase2[:,:])
                        else:
                            maskDP[:,:] = (maskBase1[:,:] | maskBase2[:,:])
                            
    # Done!
    h.close()


if __name__ == "__main__":
    main(sys.argv[1:])
    
