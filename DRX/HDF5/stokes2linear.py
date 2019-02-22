#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Given an HDF5 file with Stokes parameterss, calculate the linear polarization
products and save the results to a new file.

$Rev$
$LastChangedBy$
$LastChangedDate$
"""

import os
import sys
import h5py
import time
import numpy
import shutil
import argparse
from datetime import datetime


def _cloneStructure(input, output, level=0):
    """
    Function to recursively copy the structure of a HDF5 file created by 
    hdfWaterfall.py or drspec2hdf.py.
    """
    
    # Copy the attributes
    for key in input.attrs:
        output.attrs[key] = input.attrs[key]
        
    # Loop over the entities in the first input file for copying
    for ent in list(input):
        ## Get the entity
        entity = input.get(ent, None)
        print "%sCopying '%s'..." % (' '*level*2, ent)
        
        ## Is it a group?
        if type(entity).__name__ == 'Group':
            ### If so, add it and fill it in.
            if ent not in list(output):
                entityO = output.create_group(ent)
            _cloneStructure(entity, entityO, level=level+1)
            continue
            
        ## Is it a dataset?
        if type(entity).__name__ == 'Dataset':
            ### If so, add it and fill it in
            if ent in ('Steps', 'Delays', 'Gains', 'time', 'Saturation', 'freq'):
                entityO = output.create_dataset(ent, entity.shape, entity.dtype.descr[0][1])
                entityO[:] = entity[:]
                
            else:
                ### Don't copy the spectral data since we want to 'adjust' it
                print "%sSpectral data - skipping" % (' '*(level+1)*2,)
                continue
                
            ### Update the dataset attributes
            for key in entity.attrs:
                entityO.attrs[key] = entity.attrs[key]
                
    return True


def main(args):
    # Open the input file
    hIn = h5py.File(args.filename, mode='a')
    
    # Is it already a Stokes file and does it have enough information to continue?
    for obsName in hIn.keys():
        obs = hIn.get(obsName, None)
        tuning = obs.get('Tuning1', None)
        dataProducts = list(tuning)
        for exclude in ('freq', 'Mask', 'Saturation', 'SpectralKurtosis'):
            try:
                ind = dataProducts.index(exclude)
                del dataProducts[ind]
            except ValueError:
                pass
                
        ## Linear
        if 'XX' in dataProducts or 'YY' in dataProducts or 'XY' in dataProducts or 'YX' in dataProducts:
            raise RuntimeError("Input file '%s' contains data products %s" % (args.filename, ', '.join(dataProducts)))
            
        ## Minimum information
        pairs = 0
        if 'I' in dataProducts and 'Q' in dataProducts:
            pairs += 1
        if 'U' in dataProducts and 'V' in dataProducts:
            pairs += 1
        if pairs == 0:
            raise RuntimeError("Input file '%s' contains too few data products to form linear products" % args.filename)
            
    # Ready the output filename
    outname = os.path.basename(args.filename)
    outname = os.path.splitext(outname)[0]
    outname = "%s-linear.hdf5" % outname
    
    if os.path.exists(outname):
        if not args.force:
            yn = raw_input("WARNING: '%s' exists, overwrite? [Y/n] " % outname)
        else:
            yn = 'y'
            
        if yn not in ('n', 'N'):
            os.unlink(outname)
        else:
            raise RuntimeError("Output file '%s' already exists" % outname)
            
    hOut = h5py.File(outname, mode='w')
    
    # Copy the basic structure
    _cloneStructure(hIn, hOut)
    print " "
    
    # Loop over the observation sections
    for obsName in hIn.keys():
        obsIn = hIn.get(obsName, None)
        obsOut = hOut.get(obsName, None)
        
        # Load in the information we need
        tInt = obsIn.attrs['tInt']
        LFFT = obsIn.attrs['LFFT']
        srate = obsIn.attrs['sampleRate']
        
        print "Staring Observation #%i" % int(obsName.replace('Observation', ''))
        print "  Sample Rate: %.1f Hz" % srate
        print "  LFFT: %i" % LFFT
        print "  tInt: %.3f s" % tInt
        
        for tuning in (1, 2):
            print "    Tuning %i" % tuning
            tuningIn = obsIn.get('Tuning%i' % tuning, None)
            tuningOut = obsOut.get('Tuning%i' % tuning, None)
            
            baseMaskIn = tuningIn.get('Mask', None)
            baseMaskOut = tuningOut.get('Mask', None)
            baseSKIn = tuningIn.get('SpectralKurtosis', None)
            baseSKOut = tuningOut.get('SpectralKurtosis', None)
            
            dataProducts = list(tuningIn)
            for exclude in ('freq', 'Mask', 'Saturation', 'SpectralKurtosis'):
                try:
                    ind = dataProducts.index(exclude)
                    del dataProducts[ind]
                except ValueError:
                    pass
                    
            if 'I' in dataProducts and 'Q' in dataProducts:
                ## XX
                print "      Computing 'XX'"
                tuningOut.create_dataset('XX', tuningIn['I'].shape, dtype=tuningIn['I'].dtype.descr[0][1])
                tuningOut['XX'][:] = (tuningIn['I'][:] + tuningIn['Q'])/2.0
                
                for key in tuningIn['I'].attrs:
                    tuningOut['XX'].attrs[key] = tuningIn['I'].attrs[key]
                    
                if baseMaskIn is not None:
                    baseMaskOut.create_dataset('XX', baseMaskIn['I'].shape, dtype=baseMaskIn['I'].dtype.descr[0][1])
                    baseMaskOut['XX'][:] = baseMaskIn['I'][:] & baseMaskIn['Q'][:]
                    
                    for key in baseMaskIn['I'].attrs:
                        baseMaskOut['XX'].attrs[key] = baseMaskIn['I'].attrs[key]
                        
                if baseSKIn is not None:
                    baseSKOut.create_dataset('XX', baseSKIn['I'].shape, dtype=baseSKIn['I'].dtype.descr[0][1])
                    baseSKOut['XX'][:] = (baseSKIn['I'][:] + baseSKIn['Q'][:]) / 2.0
                    
                    for key in baseSKIn['I'].attrs:
                        baseSKOut['XX'].attrs[key] = baseSKIn['I'].attrs[key]
                        
                ## YY
                print "      Computing 'YY'"
                tuningOut.create_dataset('YY', tuningIn['I'].shape, dtype=tuningIn['I'].dtype.descr[0][1])
                tuningOut['YY'][:] = (tuningIn['I'][:] - tuningIn['Q'])/2.0
                
                for key in tuningIn['I'].attrs:
                    tuningOut['YY'].attrs[key] = tuningIn['I'].attrs[key]
                    
                if baseMaskIn is not None:
                    baseMaskOut.create_dataset('YY', baseMaskIn['I'].shape, dtype=baseMaskIn['I'].dtype.descr[0][1])
                    baseMaskOut['YY'][:] = baseMaskIn['I'][:] & baseMaskIn['Q'][:]
                    
                    for key in baseMaskIn['I'].attrs:
                        baseMaskOut['YY'].attrs[key] = baseMaskIn['I'].attrs[key]
                        
                if baseSKIn is not None:
                    baseSKOut.create_dataset('YY', baseSKIn['I'].shape, dtype=baseSKIn['I'].dtype.descr[0][1])
                    baseSKOut['YY'][:] = (baseSKIn['I'][:] + baseSKIn['Q'][:]) / 2.0
                    
                    for key in baseSKIn['I'].attrs:
                        baseSKOut['YY'].attrs[key] = baseSKIn['I'].attrs[key]
                        
            if 'U' in dataProducts and 'V' in dataProducts:
                ## XY
                print "      Computing 'XY'"
                tuningOut.create_dataset('XY', tuningIn['U'].shape, dtype=tuningIn['U'].dtype.descr[0][1])
                tuningOut['XY'][:] = (tuningIn['U'][:] + tuningIn['V'])/2.0
                
                for key in tuningIn['U'].attrs:
                    tuningOut['XY'].attrs[key] = tuningIn['U'].attrs[key]
                    
                if baseMaskIn is not None:
                    baseMaskOut.create_dataset('XY', baseMaskIn['U'].shape, dtype=baseMaskIn['U'].dtype.descr[0][1])
                    baseMaskOut['XY'][:] = baseMaskIn['U'][:] & baseMaskIn['V'][:]
                    
                    for key in baseMaskIn['U'].attrs:
                        baseMaskOut['XY'].attrs[key] = baseMaskIn['U'].attrs[key]
                        
                if baseSKIn is not None:
                    baseSKOut.create_dataset('XY', baseSKIn['U'].shape, dtype=baseSKIn['U'].dtype.descr[0][1])
                    baseSKOut['XY'][:] = (baseSKIn['U'][:] + baseSKIn['V'][:]) / 2.0
                    
                    for key in baseSKIn['U'].attrs:
                        baseSKOut['XY'].attrs[key] = baseSKIn['U'].attrs[key]
                        
                ## YX
                print "      Computing 'YX'"
                tuningOut.create_dataset('YX', tuningIn['U'].shape, dtype=tuningIn['U'].dtype.descr[0][1])
                tuningOut['YX'][:] = (tuningIn['U'][:] - tuningIn['V'])/2.0
                
                for key in tuningIn['U'].attrs:
                    tuningOut['YX'].attrs[key] = tuningIn['U'].attrs[key]
                    
                if baseMaskIn is not None:
                    baseMaskOut.create_dataset('YX', baseMaskIn['U'].shape, dtype=baseMaskIn['U'].dtype.descr[0][1])
                    baseMaskOut['YX'][:] = baseMaskIn['U'][:] & baseMaskIn['V'][:]
                    
                    for key in baseMaskIn['U'].attrs:
                        baseMaskOut['YX'].attrs[key] = baseMaskIn['U'].attrs[key]
                        
                if baseSKIn is not None:
                    baseSKOut.create_dataset('YX', baseSKIn['U'].shape, dtype=baseSKIn['U'].dtype.descr[0][1])
                    baseSKOut['YX'][:] = (baseSKIn['U'][:] + baseSKIn['V'][:]) / 2.0
                    
                    for key in baseSKIn['U'].attrs:
                        baseSKOut['YX'].attrs[key] = baseSKIn['U'].attrs[key]
    # Done!
    hIn.close()
    hOut.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='read in DRX/HDF5 waterfall file and convert the Stokes parameters into linear polarization products', 
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
    parser.add_argument('filename', type=str, 
                        help='filename to decimate')
    parser.add_argument('-f', '--force', action='store_true', 
                        help='force overwritting of existing HDF5 files')
    args = parser.parse_args()
    main(args)
    
