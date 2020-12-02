#!/usr/bin/env python

"""
Given an HDF5 file with linear polarization products, calculate the Stokes 
parameters and save the results to a new file.
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
                
        ## Stokes
        if 'I' in data_products or 'Q' in data_products or 'U' in data_products or 'V' in data_products:
            raise RuntimeError("Input file '%s' contains data products %s" % (args.filename, ', '.join(data_products)))
            
        ## Minimum information
        pairs = 0
        if 'XX' in data_products and 'YY' in data_products:
            pairs += 1
        if 'XY' in data_products and 'YX' in data_products:
            pairs += 1
        elif 'XY_real' in data_products and 'XY_imag' in data_products:
            pairs += 1
        if pairs == 0:
            raise RuntimeError("Input file '%s' contains too few data products to form Stokes parameters" % args.filename)
            
    # Ready the output filename
    outname = os.path.basename(args.filename)
    outname = os.path.splitext(outname)[0]
    outname = "%s-stokes.hdf5" % outname
    
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
            if tuningIn is None:
                continue
                
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
                    
            if 'XX' in data_products and 'YY' in data_products:
                ## I
                print("      Computing 'I'")
                tuningOut.create_dataset('I', tuningIn['XX'].shape, dtype=tuningIn['XX'].dtype.descr[0][1])
                tuningOut['I'][:] = tuningIn['XX'][:] + tuningIn['YY']
                
                for key in tuningIn['XX'].attrs:
                    tuningOut['I'].attrs[key] = tuningIn['XX'].attrs[key]
                    
                if baseMaskIn is not None:
                    baseMaskOut.create_dataset('I', baseMaskIn['XX'].shape, dtype=baseMaskIn['XX'].dtype.descr[0][1])
                    baseMaskOut['I'][:] = baseMaskIn['XX'][:] & baseMaskIn['YY'][:]
                    
                    for key in baseMaskIn['XX'].attrs:
                        baseMaskOut['I'].attrs[key] = baseMaskIn['XX'].attrs[key]
                        
                if baseSKIn is not None:
                    baseSKOut.create_dataset('I', baseSKIn['XX'].shape, dtype=baseSKIn['XX'].dtype.descr[0][1])
                    baseSKOut['I'][:] = (baseSKIn['XX'][:] + baseSKIn['YY'][:]) / 2.0
                    
                    for key in baseSKIn['XX'].attrs:
                        baseSKOut['I'].attrs[key] = baseSKIn['XX'].attrs[key]
                        
                ## Q
                print("      Computing 'Q'")
                tuningOut.create_dataset('Q', tuningIn['XX'].shape, dtype=tuningIn['XX'].dtype.descr[0][1])
                tuningOut['Q'][:] = tuningIn['XX'][:] - tuningIn['YY']
                
                for key in tuningIn['XX'].attrs:
                    tuningOut['Q'].attrs[key] = tuningIn['XX'].attrs[key]
                    
                if baseMaskIn is not None:
                    baseMaskOut.create_dataset('Q', baseMaskIn['XX'].shape, dtype=baseMaskIn['XX'].dtype.descr[0][1])
                    baseMaskOut['Q'][:] = baseMaskIn['XX'][:] & baseMaskIn['YY'][:]
                    
                    for key in baseMaskIn['XX'].attrs:
                        baseMaskOut['Q'].attrs[key] = baseMaskIn['XX'].attrs[key]
                        
                if baseSKIn is not None:
                    baseSKOut.create_dataset('Q', baseSKIn['XX'].shape, dtype=baseSKIn['XX'].dtype.descr[0][1])
                    baseSKOut['Q'][:] = (baseSKIn['XX'][:] + baseSKIn['YY'][:]) / 2.0
                    
                    for key in baseSKIn['XX'].attrs:
                        baseSKOut['Q'].attrs[key] = baseSKIn['XX'].attrs[key]
                        
            if 'XY' in data_products and 'YX' in data_products:
                ## U
                print("      Computing 'U'")
                tuningOut.create_dataset('U', tuningIn['XY'].shape, dtype=tuningIn['XY'].dtype.descr[0][1])
                tuningOut['U'][:] = tuningIn['XY'][:] + tuningIn['YX']
                
                for key in tuningIn['XY'].attrs:
                    tuningOut['U'].attrs[key] = tuningIn['XY'].attrs[key]
                    
                if baseMaskIn is not None:
                    baseMaskOut.create_dataset('U', baseMaskIn['XY'].shape, dtype=baseMaskIn['XY'].dtype.descr[0][1])
                    baseMaskOut['U'][:] = baseMaskIn['XY'][:] & baseMaskIn['YX'][:]
                    
                    for key in baseMaskIn['XY'].attrs:
                        baseMaskOut['U'].attrs[key] = baseMaskIn['XY'].attrs[key]
                        
                if baseSKIn is not None:
                    baseSKOut.create_dataset('U', baseSKIn['XY'].shape, dtype=baseSKIn['XY'].dtype.descr[0][1])
                    baseSKOut['U'][:] = (baseSKIn['XY'][:] + baseSKIn['YX'][:]) / 2.0
                    
                    for key in baseSKIn['XY'].attrs:
                        baseSKOut['U'].attrs[key] = baseSKIn['XY'].attrs[key]
                        
                ## V
                print("      Computing 'V'")
                tuningOut.create_dataset('V', tuningIn['XY'].shape, dtype=tuningIn['XY'].dtype.descr[0][1])
                tuningOut['V'][:] = tuningIn['XY'][:] - tuningIn['YX']
                
                for key in tuningIn['XY'].attrs:
                    tuningOut['V'].attrs[key] = tuningIn['XY'].attrs[key]
                    
                if baseMaskIn is not None:
                    baseMaskOut.create_dataset('V', baseMaskIn['XY'].shape, dtype=baseMaskIn['XY'].dtype.descr[0][1])
                    baseMaskOut['V'][:] = baseMaskIn['XY'][:] & baseMaskIn['YX'][:]
                    
                    for key in baseMaskIn['XY'].attrs:
                        baseMaskOut['V'].attrs[key] = baseMaskIn['XY'].attrs[key]
                        
                if baseSKIn is not None:
                    baseSKOut.create_dataset('V', baseSKIn['XY'].shape, dtype=baseSKIn['XY'].dtype.descr[0][1])
                    baseSKOut['V'][:] = (baseSKIn['XY'][:] + baseSKIn['YX'][:]) / 2.0
                    
                    for key in baseSKIn['XY'].attrs:
                        baseSKOut['V'].attrs[key] = baseSKIn['XY'].attrs[key]
                        
            elif 'XY_real' in data_products and 'XY_imag' in data_products:
                XY = tuningIn['XY_real'][:] + 1j*tuningIn['XY_imag'][:]
                
                ## U
                print("      Computing 'U'")
                tuningOut.create_dataset('U', tuningIn['XY_real'].shape, dtype=tuningIn['XY_real'].dtype.descr[0][1])
                tuningOut['U'][:] = 2*XY.real
                
                for key in tuningIn['XY_real'].attrs:
                    tuningOut['U'].attrs[key] = tuningIn['XY_real'].attrs[key]
                    
                if baseMaskIn is not None:
                    baseMaskOut.create_dataset('U', baseMaskIn['XY_real'].shape, dtype=baseMaskIn['XY_real'].dtype.descr[0][1])
                    baseMaskOut['U'][:] = baseMaskIn['XY_real'][:] & baseMaskIn['XY_imag'][:]
                    
                    for key in baseMaskIn['XY_real'].attrs:
                        baseMaskOut['U'].attrs[key] = baseMaskIn['XY_real'].attrs[key]
                        
                if baseSKIn is not None:
                    baseSKOut.create_dataset('U', baseSKIn['XY_real'].shape, dtype=baseSKIn['XY_real'].dtype.descr[0][1])
                    baseSKOut['U'][:] = (baseSKIn['XY_real'][:] + baseSKIn['XY_imag'][:]) / 2.0
                    
                    for key in baseSKIn['XY_real'].attrs:
                        baseSKOut['U'].attrs[key] = baseSKIn['XY_real'].attrs[key]
                        
                ## V
                print("      Computing 'V'")
                tuningOut.create_dataset('V', tuningIn['XY_real'].shape, dtype=tuningIn['XY_real'].dtype.descr[0][1])
                tuningOut['V'][:] = 2*XY.imag
                
                for key in tuningIn['XY_real'].attrs:
                    tuningOut['V'].attrs[key] = tuningIn['XY_real'].attrs[key]
                    
                if baseMaskIn is not None:
                    baseMaskOut.create_dataset('V', baseMaskIn['XY_real'].shape, dtype=baseMaskIn['XY_real'].dtype.descr[0][1])
                    baseMaskOut['V'][:] = baseMaskIn['XY_real'][:] & baseMaskIn['XY_imag'][:]
                    
                    for key in baseMaskIn['XY_real'].attrs:
                        baseMaskOut['V'].attrs[key] = baseMaskIn['XY_real'].attrs[key]
                        
                if baseSKIn is not None:
                    baseSKOut.create_dataset('V', baseSKIn['XY_real'].shape, dtype=baseSKIn['XY_real'].dtype.descr[0][1])
                    baseSKOut['V'][:] = (baseSKIn['XY_real'][:] + baseSKIn['XY_imag'][:]) / 2.0
                    
                    for key in baseSKIn['XY_real'].attrs:
                        baseSKOut['V'].attrs[key] = baseSKIn['XY_real'].attrs[key]
                        
    # Done!
    hIn.close()
    hOut.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='read in DRX/HDF5 waterfall file and convert the linear polarization products to Stokes parameters', 
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
    parser.add_argument('filename', type=str, 
                        help='filename to decimate')
    parser.add_argument('-f', '--force', action='store_true', 
                        help='force overwritting of existing HDF5 files')
    args = parser.parse_args()
    main(args)
    
