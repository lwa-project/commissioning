#!/usr/bin/env python

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
"""

# Python3 compatiability
from __future__ import print_function, division
import sys
if sys.version_info > (3,):
    xrange = range
    
import os
import sys
import h5py
import numpy
import argparse
from datetime import datetime

from lsl.statistics import robust, kurtosis
from lsl.misc import parser as aph


def main(args):
    # Open the file and loop over the observation sections
    h = h5py.File(args.filename, mode='a' if (not args.no_update) else 'r')
    for obsName in h.keys():
        obs = h.get(obsName, None)
        
        # Load in the information we need to calculate the pseudo-spectral kurtosis (pSK)
        tInt = obs.attrs['tInt']
        LFFT = obs.attrs['LFFT']
        srate = obs.attrs['sampleRate']
        
        skN = int(tInt*srate / LFFT)
        chunkSize = int(round(args.duration/tInt))
        
        print("Staring Observation #%i" % int(obsName.replace('Observation', '')))
        print("  Sample Rate: %.1f Hz" % srate)
        print("  LFFT: %i" % LFFT)
        print("  tInt: %.3f s" % tInt)
        print("  Integrations per spectrum: %i" % skN)
        
        time = obs.get('time', None)[:]
        tuning1 = obs.get('Tuning1', None)
        tuning2 = obs.get('Tuning2', None)
        
        # Get all of the avaliable data products
        data_products = list(tuning1)
        for toRemove in ('Mask', 'Saturation', 'SpectralKurtosis', 'freq'):
            try:
                del data_products[data_products.index(toRemove)]
            except ValueError:
                pass
                
        # What's this?  Stokes I is the sum of XX and YY so there are 
        # effectively twice as many FFTs per spectrum relative to XX and
        # YY
        nAdjust = {'XX': 1, 'YY': 1, 'I': 2, 'RR': 1, 'LL': 1}
        
        # Loop over data products - primary
        for dp in ('XX', 'YY', 'I', 'RR', 'LL'):
            for t,tuning in enumerate((tuning1, tuning2)):
                if tuning is None:
                    continue
                    
                data = tuning.get(dp, None)
                if data is None:
                    continue
                    
                # Loop over chunks for computing pSK
                sk = numpy.zeros_like(data)
                
                section = numpy.empty((chunkSize,data.shape[1]), dtype=data.dtype)
                for i in xrange(time.size//chunkSize+1):
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
                        sk[start:stop,j] = kurtosis.spectral_power(section[:,j], N=skN*nAdjust[dp])
                        
                # Report
                print("  => %s-%i SK Mean: %.3f" % (dp, t+1, numpy.mean(sk)))
                print("     %s-%i SK Std. Dev.: %.3f" % (dp, t+1, numpy.std(sk)))
                
                # Save the pSK information to the HDF5 file if we need to
                if (not args.no_update):
                    h.attrs['FileLastUpdated'] = datetime.utcnow().strftime("UTC %Y/%m/%d %H:%M:%S")
                    
                    if args.generate_mask:
                        ## Calculate the expected pSK limits for the threshold
                        kLow, kHigh = kurtosis.get_limits(args.threshold, chunkSize, N=skN*nAdjust[dp])
                        
                        ## Generate the mask arrays
                        maskTable = numpy.where( (sk < kLow) | (sk > kHigh), True, False )
                        
                        ## Report
                        print("     => %s-%i Mask Fraction: %.1f%%" % (dp, t+1, 100.0*maskTable.sum()/maskTable.size,))
                        
                        ## Pull out the Mask group from the HDF5 file
                        mask = tuning.get('Mask', None)
                        if mask is None:
                            mask = tuning.create_group('Mask')
                            for p in data_products:
                                mask.create_dataset(p, tuning[p].shape, 'bool')
                        maskDP = mask.get(dp, None)
                        if args.merge:
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
                        
        if args.generate_mask and args.fill:	
            # Loop over data products - secondary
            for dp in data_products:
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
                        if args.merge:
                            maskDP[:,:] |= maskBase[:,:]
                        else:
                            maskDP[:,:] = maskBase[:,:]
                            
                    elif dp in ('XY_real', 'XY_imag'):
                        ## Case 2 - Flag XY_real and XY_imag off XX and YY
                        
                        ## Pull out the Mask group from the HDF5 file
                        mask = tuning.get('Mask', None)
                        maskDP = mask.get(dp, None)
                        maskBase1 = mask.get('XX', None)
                        maskBase2 = mask.get('YY', None)
                        if args.merge:
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
                        if args.merge:
                            maskDP[:,:] |= (maskBase1[:,:] | maskBase2[:,:])
                        else:
                            maskDP[:,:] = (maskBase1[:,:] | maskBase2[:,:])
                            
    # Done!
    h.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='read in DRX/HDF5 waterfall file and calculate the pseudo-spectral kurtosis (pSK)', 
        epilog='NOTE:  The -g/--generate-mask flag overrides the -n/--no-update flag.', 
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
    parser.add_argument('filename', type=str, 
                        help='filename to process')
    parser.add_argument('-d', '--duration', type=aph.positive_float, default=10.0, 
                        help='pSK update interval in seconds')
    parser.add_argument('-n', '--no-update', action='store_true', 
                        help='do not add the pSK information to the HDF5 file')
    parser.add_argument('-g', '--generate-mask', action='store_true', 
                        help='add pSK information to the file as a mask rather than pSK values')
    parser.add_argument('-t', '--threshold', type=aph.positive_float, default=4.0, 
                        help='masking threshold, in sigma, if the -g/--generate-mask flag is set')
    parser.add_argument('-f', '--fill', action='store_true', 
                        help='fill in the secondary polarizations, i.e., XY_real, XY_imag; Q, U, V; RL, LR, using the mask from the primaries: XX, YY; I; RR, LL')
    parser.add_argument('-m', '--merge', action='store_true', 
                        help='merge the new mask table with the existing one, if it exists')
    args = parser.parse_args()
    main(args)
    
