#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Export select stands from a TBW file on a DRSU to HDF5.

Usage:
./tbw2hdfDRSU.py <drsu_device> <drsu_TBW_tag> <stand_ID> [<stand_ID> [...]]

$Rev$
$LastChangedBy$
$LastChangedDate$
"""

import os
import sys
import h5py
import ephem
import numpy

from lsl.common.stations import lwa1
from lsl.reader import drsu, tbw, errors
from lsl.astro import unix_to_utcjd, DJD_OFFSET


def main(args):
	device = args[0]
	filename = args[1]
	stands = [int(i) for i in args[2:]]

	antennas = lwa1.getAntennas()
	
	# Build the TBW file
	tbwFile = drsu.getFileByName(device, filename)

	tbwFile.open()
	nFrames = tbwFile.size / tbw.FrameSize
	dataBits = tbw.getDataBits(tbwFile.fh)
	# The number of ant/pols in the file is hard coded because I cannot figure out 
	# a way to get this number in a systematic fashion
	antpols = len(antennas)
	if dataBits == 12:
		nSamples = 400
	else:
		nSamples = 1200

	# Read in the first frame and get the date/time of the first sample 
	# of the frame.  This is needed to get the list of stands.
	junkFrame = tbw.readFrame(tbwFile.fh)
	tbwFile.seek(-tbw.FrameSize, 1)
	beginTime = junkFrame.getTime()
	beginDate = ephem.Date(unix_to_utcjd(junkFrame.getTime()) - DJD_OFFSET)
	
	# Figure out which digitizers to keep
	toKeep = []
	for a in antennas:
		if a.stand.id in stands:
			toKeep.append( a.digitizer )

	# File summary
	print "Filename: %s" % filename
	print "Date of First Frame: %s" % str(beginDate)
	print "Ant/Pols: %i" % antpols
	print "Sample Length: %i-bit" % dataBits
	print "Frames: %i" % nFrames
	print "==="
	print "Keeping Stands:"
	for a in toKeep:
		print " Stand #%3i, pol %i (digitizer %3i)" % (antennas[a-1].stand.id, antennas[a-1].pol, antennas[a-1].digitizer)

	# Skip over any non-TBW frames at the beginning of the file
	i = 0
	junkFrame = tbw.readFrame(tbwFile.fh)
	while not junkFrame.header.isTBW():
		junkFrame = tbw.readFrame(tbwFile.fh)
		i += 1
	tbwFile.seek(-tbw.FrameSize, 1)
	print "Skipped %i non-TBW frames at the beginning of the file" % i
	
	# Create the HDF5 file
	outname = "%s_TBW.hdf5" % filename
	f = h5py.File(outname, 'w')
	f.attrs['filename'] = filename
	f.attrs['mode'] = 'TBW'
	f.attrs['station'] = 'LWA-1'
	f.attrs['dataBits'] = dataBits
	f.attrs['startTime'] = beginTime
	f.attrs['startTime_units'] = 's'
	f.attrs['startTime_sys'] = 'unix'
	f.attrs['sampleRate'] = 196e6
	f.attrs['sampleRate_units'] = 'Hz'
	
	## Create the digitzer to dataset lookup table and the 
	standLookup = {}
	standData = []
	i = 0
	for a in toKeep:
		if a % 2 == 0:
			continue
		
		s = antennas[a-1].stand.id
		standLookup[a] = i
		
		temp = f.create_group('Stand%03i' % s)
		### Combined status code
		temp.attrs['statusCode'] = antennas[a-1].getStatus()
		
		### Antenna number
		temp.attrs['antennaID'] = antennas[a-1].id
		
		### Cable information
		temp.attrs['cableID'] = antennas[a-1].cable.id
		temp.attrs['cableLength'] = antennas[a-1].cable.length
		temp.attrs['cableLength_units'] = 'm'
		
		### Stand location information
		temp.attrs['posX'] = antennas[a-1].stand.x
		temp.attrs['posX_units'] = 'm'
		temp.attrs['posY'] = antennas[a-1].stand.y
		temp.attrs['posY_units'] = 'm'
		temp.attrs['posZ'] = antennas[a-1].stand.z
		temp.attrs['posZ_units'] = 'm'
		
		### Time series data sets
		temp.attrs['axis0'] = 'time'
		xpol = temp.create_dataset('X', (12000000,), 'i2', chunks=True)
		ypol = temp.create_dataset('Y', (12000000,), 'i2', chunks=True)
		
		standData.append( (xpol, ypol, temp) )
		i += 1
	
	# Go!
	while True:
		# Read in the next frame and anticipate any problems that could occur
		try:
			cFrame = tbw.readFrame(tbwFile.fh)
		except errors.eofError:
			break
		except errors.syncError:
			print "WARNING: Mark 5C sync error on frame #%i" % (int(tbwFile.tell())/tbw.FrameSize-1)
			continue
		if not cFrame.header.isTBW():
			continue
		
		# Get the DP "stand" ID and the digitizer number
		stand = cFrame.header.parseID()
		aStand = 2*(stand-1)
		digitizer = aStand + 1
		
		# If we don't need it, skip it
		if digitizer not in toKeep:
			continue
		
		# Actually load the data.
		## Frame count
		count = cFrame.header.frameCount - 1
		## Which data set
		dataset = standLookup[digitizer]
		## Load
		standData[dataset][0][count*nSamples:(count+1)*nSamples] = cFrame.data.xy[0,:]
		standData[dataset][1][count*nSamples:(count+1)*nSamples] = cFrame.data.xy[1,:]
		
	tbwFile.close()
	f.close()


if __name__ == "__main__":
	main(sys.argv[1:])
	