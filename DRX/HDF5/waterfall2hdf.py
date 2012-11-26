#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Convert the output NPZ of drxWaterfall/drxWaterfallDRSU into a HDF5 file.

$Rev$
$LastChangedBy$
$LastChangedDate$
"""

import os
import sys
import h5py
import numpy


def main(args):
	for filename in args:
		dataDict = numpy.load(filename)
		
		freq  = dataDict['freq']
		freq1 = dataDict['freq1']
		freq2 = dataDict['freq2']
		
		tInt = dataDict['tInt']
		times = dataDict['times']
		sRate = dataDict['srate']
		standMapper = dataDict['standMapper']
		
		spec = dataDict['spec']
		try:
			mask = dataDict['mask'].astype(numpy.bool)
		except KeyError:
			mask = None
		
		hdfFilename = filename.replace('-waterfall.npz', '-waterfall.hdf5')
		f = h5py.File(hdfFilename, 'w')
		f.attrs['Beam'] = standMapper[0]/4 + 1
		f.attrs['tInt'] = tInt
		f.attrs['tInt_Units'] = 's'
		f.attrs['sampleRate'] = sRate
		f.attrs['sampleRate_Units'] = 'samples/s'
		f.attrs['RBW'] = freq[1]-freq[0]
		f.attrs['RBW_Units'] = 'Hz'
		f['time'] = times.astype(numpy.float64)
		
		tuning1 = f.create_group('/Tuning1')
		tuning1['freq'] = freq1.astype(numpy.float64)
		tuning1['freq'].attrs['Units'] = 'Hz'
		tuning1['XX'] = numpy.squeeze(spec[:,0,:].astype(numpy.float32))
		tuning1['XX'].chunks = True
		tuning1['XX'].attrs['axis0'] = 'time'
		tuning1['XX'].attrs['axis1'] = 'frequency'
		tuning1['YY'] = numpy.squeeze(spec[:,1,:].astype(numpy.float32))
		tuning1['YY'].chunks = True
		tuning1['YY'].attrs['axis0'] = 'time'
		tuning1['YY'].attrs['axis1'] = 'frequency'
		
		tuning2 = f.create_group('/Tuning2')
		tuning2['freq'] = freq2.astype(numpy.float64)
		tuning2['freq'].attrs['Units'] = 'Hz'
		tuning2['XX'] = numpy.squeeze(spec[:,2,:].astype(numpy.float32))
		tuning2['XX'].chunks = True
		tuning2['XX'].attrs['axis0'] = 'time'
		tuning2['XX'].attrs['axis1'] = 'frequency'
		tuning2['YY'] = numpy.squeeze(spec[:,3,:].astype(numpy.float32))
		tuning2['YY'].chunks = True
		tuning2['YY'].attrs['axis0'] = 'time'
		tuning2['YY'].attrs['axis1'] = 'frequency'
		
		if mask is not None:
			mask1 = tuning1.create_group('Mask')
			mask1X = mask1.create_dataset('XX', tuning1['XX'].shape, 'bool', chunks=True)
			mask1X = mask[:,0,:]
			mask1Y = mask1.create_dataset('YY', tuning1['YY'].shape, 'bool', chunks=True)
			mask1Y = mask[:,1,:]
			
			mask2 = tuning2.create_group('Mask')
			mask2X = mask2.create_dataset('XX', tuning2['XX'].shape, 'bool', chunks=True)
			mask2X = mask[:,2,:]
			mask2Y = mask2.create_dataset('YY', tuning2['YY'].shape, 'bool', chunks=True)
			mask2Y = mask[:,3,:]
		
		f.close()
		
		
if __name__ == "__main__":
	main(sys.argv[1:])
	