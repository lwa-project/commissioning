#!/usr/bin/env python

"""
Given an HDF5 file with Stokes parameterss, calculate the linear polarization
products and save the results to a new file.
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
        print("%sCopying '%s'..." % (' '*level*2, ent))
        
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
                entityO = output.create_dataset(ent, entity.shape, entity.dtype)
                entityO[:] = entity[:]
                
            else:
                ### Don't copy the spectral data since we want to 'adjust' it
                print("%sSpectral data - skipping" % (' '*(level+1)*2,))
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
        data_products = list(tuning)
        for exclude in ('freq', 'Mask', 'Saturation', 'SpectralKurtosis'):
            try:
                ind = data_products.index(exclude)
                del data_products[ind]
            except ValueError:
                pass
                
        ## Linear
        if 'XX' in data_products or 'YY' in data_products or 'XY_real' in data_products or 'XY_imag' in data_products:
            raise RuntimeError("Input file '%s' contains data products %s" % (args.filename, ', '.join(data_products)))
            
        ## Minimum information
        pairs = 0
        if 'I' in data_products and 'Q' in data_products:
            pairs += 1
        if 'U' in data_products and 'V' in data_products:
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
    print(" ")
    
    # Loop over the observation sections
    for obsName in hIn.keys():
        obsIn = hIn.get(obsName, None)
        obsOut = hOut.get(obsName, None)
        
        # Load in the information we need
        tInt = obsIn.attrs['tInt']
        LFFT = obsIn.attrs['LFFT']
        srate = obsIn.attrs['sampleRate']
        
        print("Staring Observation #%i" % int(obsName.replace('Observation', '')))
        print("  Sample Rate: %.1f Hz" % srate)
        print("  LFFT: %i" % LFFT)
        print("  tInt: %.3f s" % tInt)
        
        for tuning in (1, 2):
            print("    Tuning %i" % tuning)
            tuningIn = obsIn.get('Tuning%i' % tuning, None)
            tuningOut = obsOut.get('Tuning%i' % tuning, None)
            
            baseMaskIn = tuningIn.get('Mask', None)
            baseMaskOut = tuningOut.get('Mask', None)
            baseSKIn = tuningIn.get('SpectralKurtosis', None)
            baseSKOut = tuningOut.get('SpectralKurtosis', None)
            
            data_products = list(tuningIn)
            for exclude in ('freq', 'Mask', 'Saturation', 'SpectralKurtosis'):
                try:
                    ind = data_products.index(exclude)
                    del data_products[ind]
                except ValueError:
                    pass
                    
            if 'I' in data_products and 'Q' in data_products:
                ## XX
                print("      Computing 'XX'")
                tuningOut.create_dataset('XX', tuningIn['I'].shape, dtype=tuningIn['I'].dtype.descr[0][1])
                tuningOut['XX'][:] = (tuningIn['I'][:] + tuningIn['Q'][:])/2.0
                
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
                print("      Computing 'YY'")
                tuningOut.create_dataset('YY', tuningIn['I'].shape, dtype=tuningIn['I'].dtype.descr[0][1])
                tuningOut['YY'][:] = (tuningIn['I'][:] - tuningIn['Q'][:])/2.0
                
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
                        
            if 'U' in data_products and 'V' in data_products:
                ## XY_real
                print("      Computing 'XY_real'")
                tuningOut.create_dataset('XY_real', tuningIn['U'].shape, dtype=tuningIn['U'].dtype.descr[0][1])
                tuningOut['XY_real'][:] = tuningIn['U'][:]/2.0
                
                for key in tuningIn['U'].attrs:
                    tuningOut['XY_real'].attrs[key] = tuningIn['U'].attrs[key]
                    
                if baseMaskIn is not None:
                    baseMaskOut.create_dataset('XY_real', baseMaskIn['U'].shape, dtype=baseMaskIn['U'].dtype.descr[0][1])
                    baseMaskOut['XY_real'][:] = baseMaskIn['U'][:]
                    
                    for key in baseMaskIn['U'].attrs:
                        baseMaskOut['XY_real'].attrs[key] = baseMaskIn['U'].attrs[key]
                        
                if baseSKIn is not None:
                    baseSKOut.create_dataset('XY_real', baseSKIn['U'].shape, dtype=baseSKIn['U'].dtype.descr[0][1])
                    baseSKOut['XY_real'][:] = baseSKIn['U'][:]
                    
                    for key in baseSKIn['U'].attrs:
                        baseSKOut['XY_real'].attrs[key] = baseSKIn['U'].attrs[key]
                        
                ## XY_imag
                print("      Computing 'XY_imag'")
                tuningOut.create_dataset('XY_imag', tuningIn['V'].shape, dtype=tuningIn['V'].dtype.descr[0][1])
                tuningOut['XY_imag'][:] = tuningIn['V'][:]/2.0
                
                for key in tuningIn['V'].attrs:
                    tuningOut['XY_imag'].attrs[key] = tuningIn['V'].attrs[key]
                    
                if baseMaskIn is not None:
                    baseMaskOut.create_dataset('XY_imag', baseMaskIn['V'].shape, dtype=baseMaskIn['V'].dtype.descr[0][1])
                    baseMaskOut['XY_imag'][:] = baseMaskIn['V'][:]
                    
                    for key in baseMaskIn['V'].attrs:
                        baseMaskOut['XY_imag'].attrs[key] = baseMaskIn['V'].attrs[key]
                        
                if baseSKIn is not None:
                    baseSKOut.create_dataset('XY_imag', baseSKIn['V'].shape, dtype=baseSKIn['V'].dtype.descr[0][1])
                    baseSKOut['XY_imag'][:] = baseSKIn['V'][:]
                    
                    for key in baseSKIn['V'].attrs:
                        baseSKOut['XY_imag'].attrs[key] = baseSKIn['V'].attrs[key]
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
    
