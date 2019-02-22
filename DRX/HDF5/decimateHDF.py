#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Given an HDF5 file, decimate the data contained in it in both time and 
frequency, and save the results to a new file.

$Rev$
$LastChangedBy$
$LastChangedDate$
"""


import os
import sys
import h5py
import numpy
import argparse

from lsl.misc import parser as aph


def _fillHDF(input, output, tDecimation=1, sDecimation=1, level=0):
    """
    Function to recursively copy the structure of a HDF5 file created by 
    hdfWaterfall.py or drspec2hdf.py.
    """
    
    # Copy the attributes
    for key in input.attrs:
        if key == 'tInt':
            value = input.attrs[key]*tDecimation
        elif key == 'nChan':
            value = input.attrs[key]/sDecimation
        elif key == 'RBW':
            value = input.attrs[key]*sDecimation
        else:
            value = input.attrs[key]
        output.attrs[key] = value
        
    # Loop over the entities in the first input file for copying
    for ent in list(input):
        ## Get the entity
        entity = input.get(ent, None)
        print "%sWorking on '%s'..." % (' '*level*2, ent)
        
        ## Is it a group?
        if type(entity).__name__ == 'Group':
            ### If so, add it and fill it in.
            if ent not in list(output):
                entityO = output.create_group(ent)
            _fillHDF(entity, entityO, tDecimation=tDecimation, sDecimation=sDecimation, level=level+1)
            continue
            
        ## Is it a dataset?
        if type(entity).__name__ == 'Dataset':
            ### If so, add it and fill it in
            if ent in ('Steps', 'Delays', 'Gains'):
                entity0 = output.create_dataset(ent, entity.shape, entity.dtype.descr[0][1])
                entity0[:] = entity[:]
                
            elif ent == 'time':
                newShape = (entity.shape[0]/tDecimation,)
                entityO = output.create_dataset(ent, newShape, entity.dtype.descr[0][1])
                for i in xrange(newShape[0]):
                    data = entity[tDecimation*i:tDecimation*(i+1)]
                    entityO[i] = data[0]
                    
            elif ent == 'Saturation':
                newShape = (entity.shape[0]/tDecimation, entity.shape[1])
                entityO = output.create_dataset(ent, newShape, entity.dtype.descr[0][1])
                for i in xrange(newShape[0]):
                    data = entity[tDecimation*i:tDecimation*(i+1),:]
                    entityO[i,:] = data.sum(axis=0)
                    
            elif ent == 'freq':
                newShape = (entity.shape[0]/sDecimation,)
                entityO = output.create_dataset(ent, newShape, entity.dtype.descr[0][1])
                for i in xrange(newShape[0]):
                    data = entity[sDecimation*i:sDecimation*(i+1)]
                    entityO[i] = data.mean()
                    
            else:
                newShape = (entity.shape[0]/tDecimation, entity.shape[1]/sDecimation)
                entityO = output.create_dataset(ent, newShape, entity.dtype.descr[0][1])
                for i in xrange(newShape[0]):
                    data = entity[tDecimation*i:tDecimation*(i+1),:newShape[1]*sDecimation]
                    data = data.mean(axis=0)
                    data.shape = (data.size/sDecimation, sDecimation)
                    data = data.mean(axis=1)
                    if data.dtype != entity.dtype:
                        data = data.astype(entity.dtype)
                    entityO[i,:] = data
                    
            ### Update the dataset attributes
            for key in entity.attrs:
                entityO.attrs[key] = entity.attrs[key]
                
    return True


def main(args):
    for filename in args.filename:
        outname = os.path.basename(filename)
        outname = os.path.splitext(outname)[0]
        outname = '%s-decim.hdf5' % outname
        
        if os.path.exists(outname):
            if not args.force:
                yn = raw_input("WARNING: '%s' exists, overwrite? [Y/n] " % outname)
            else:
                yn = 'y'
                
            if yn not in ('n', 'N'):
                os.unlink(outname)
            else:
                print "WARNING: output file '%s' already exists, skipping" % outname
                continue
                
        hIn  = h5py.File(filename, mode='r')
        hOut = h5py.File(outname, mode='w')
        
        _fillHDF(hIn, hOut, tDecimation=args.time_decimation, sDecimation=args.spec_decimation)
        
        hIn.close()
        hOut.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='read in DRX/HDF5 waterfall file and decimate the file as requested', 
        epilog='NOTE:  This scripts decimates even if the number of times steps or frquency channels is not an intger multiple of the decimation factor.  This can lead to data loss at the end of observations and at the higher channel numbers.', 
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
    parser.add_argument('time_decimation', type=aph.positive_int, 
                        help='temporal decimation factor')
    parser.add_argument('filename', type=str, nargs='+', 
                        help='filename to decimate')
    parser.add_argument('-s', '--spec-decimation', type=aph.positive_int, default=1, 
                        help='apply a decimation to the spectral dimension')
    parser.add_argument('-f', '--force', action='store_true', 
                        help='force overwritting of existing HDF5 files')
    args = parser.parse_args()
    main(args)
    
