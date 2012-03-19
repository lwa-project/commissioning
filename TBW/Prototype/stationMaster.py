#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Given a TBW file, plot the time averaged spectra for each digitizer input.  Save 
the data for later review with smGUI as an NPZ file.

$Rev$
$LastChangedBy$
$LastChangedDate$
"""

import os
import sys
import math
import numpy
import ephem
import getopt

from lsl.common import stations
from lsl.reader import tbw
from lsl.reader import errors
from lsl.correlator import fx as fxc
from lsl.astro import unix_to_utcjd, DJD_OFFSET
from lsl.common.progress import ProgressBar

import matplotlib.pyplot as plt


def usage(exitCode=None):
	print """stationMaster.py - Read in TBW files and create a collection of 
time-averaged spectra.

Usage: stationMaster.py [OPTIONS] file

Options:
-h, --help                  Display this help information
-m, --metadata              Name of SSMIF file to use for mappings
-f, --force                 Remake the NPZ file, even if it exists
-t, --bartlett              Apply a Bartlett window to the data
-b, --blackman              Apply a Blackman window to the data
-n, --hanning               Apply a Hanning window to the data
-q, --quiet                 Run tbwSpectra in silent mode
-l, --fft-length            Set FFT length (default = 4096)
"""

	if exitCode is not None:
		sys.exit(exitCode)
	else:
		return True


def parseOptions(args):
	config = {}
	# Command line flags - default values
	config['SSMIF'] = ''
	config['force'] = False
	config['LFFT'] = 4096
	config['maxFrames'] = 30000*10
	config['window'] = fxc.noWindow
	config['applyGain'] = True
	config['verbose'] = True
	config['args'] = []

	# Read in and process the command line flags
	try:
		opts, args = getopt.getopt(args, "hm:fqtbnl:", ["help", "metadata=", "force", "quiet", "bartlett", "blackman", "hanning", "fft-length="])
	except getopt.GetoptError, err:
		# Print help information and exit:
		print str(err) # will print something like "option -a not recognized"
		usage(exitCode=2)
	
	# Work through opts
	for opt, value in opts:
		if opt in ('-h', '--help'):
			usage(exitCode=0)
		elif opt in ('-m', '--metadata'):
			config['SSMIF'] = value
		elif opt in ('-f', '--force'):
			config['force'] = True
		elif opt in ('-q', '--quiet'):
			config['verbose'] = False
		elif opt in ('-t', '--bartlett'):
			config['window'] = numpy.bartlett
		elif opt in ('-b', '--blackman'):
			config['window'] = numpy.blackman
		elif opt in ('-n', '--hanning'):
			config['window'] = numpy.hanning
		elif opt in ('-l', '--fft-length'):
			config['LFFT'] = int(value)
		else:
			assert False
	
	# Add in arguments
	config['args'] = args

	# Return configuration
	return config

def main(args):
	# Parse command line options
	config = parseOptions(args)
	
	# Set the station
	if config['SSMIF'] != '':
		station = stations.parseSSMIF(config['SSMIF'])
		ssmifContents = open(config['SSMIF']).readlines()
	else:
		station = stations.lwa2
		ssmifContents = ''
	antennas = []
	for a in station.getAntennas():
		if a.digitizer != 0:
			antennas.append(a)

	# Length of the FFT
	LFFT = config['LFFT']

	# Make sure that the file chunk size contains is an integer multiple
	# of the FFT length so that no data gets dropped
	maxFrames = int(config['maxFrames']/float(LFFT))*LFFT
	# It seems like that would be a good idea, however...  TBW data comes one
	# capture at a time so doing something like this actually truncates data 
	# from the last set of stands for the first integration.  So, we really 
	# should stick with
	maxFrames = config['maxFrames']

	fh = open(config['args'][0], "rb")
	nFrames = os.path.getsize(config['args'][0]) / tbw.FrameSize
	dataBits = tbw.getDataBits(fh)
	# The number of ant/pols in the file is hard coded because I cannot figure out 
	# a way to get this number in a systematic fashion
	antpols = len(antennas)
	nChunks = int(math.ceil(1.0*nFrames/maxFrames))
	if dataBits == 12:
		nSamples = 400
	else:
		nSamples = 1200

	# Read in the first frame and get the date/time of the first sample 
	# of the frame.  This is needed to get the list of stands.
	junkFrame = tbw.readFrame(fh)
	fh.seek(0)
	beginTime = junkFrame.getTime()
	beginDate = ephem.Date(unix_to_utcjd(junkFrame.getTime()) - DJD_OFFSET)

	# File summary
	print "Filename: %s" % config['args'][0]
	print "Date of First Frame: %s" % str(beginDate)
	print "Ant/Pols: %i" % antpols
	print "Sample Length: %i-bit" % dataBits
	print "Frames: %i" % nFrames
	print "Chunks: %i" % nChunks
	print "==="

	nChunks = 1

	# Skip over any non-TBW frames at the beginning of the file
	i = 0
	junkFrame = tbw.readFrame(fh)
	while not junkFrame.header.isTBW():
		junkFrame = tbw.readFrame(fh)
		i += 1
	fh.seek(-tbw.FrameSize, 1)
	print "Skipped %i non-TBW frames at the beginning of the file" % i

	base, ext = os.path.splitext(config['args'][0])
	if (not os.path.exists("%s.npz" % base)) or config['force']:
		# Master loop over all of the file chunks
		masterSpectra = numpy.zeros((nChunks, antpols, LFFT-1))
		for i in range(nChunks):
			# Find out how many frames remain in the file.  If this number is larger
			# than the maximum of frames we can work with at a time (maxFrames),
			# only deal with that chunk
			framesRemaining = nFrames - i*maxFrames
			if framesRemaining > maxFrames:
				framesWork = maxFrames
			else:
				framesWork = framesRemaining
			print "Working on chunk %i, %i frames remaining" % ((i+1), framesRemaining)

			data = numpy.zeros((antpols, 2*30000*10*nSamples/antpols), dtype=numpy.int16)
			# If there are fewer frames than we need to fill an FFT, skip this chunk
			if data.shape[1] < 2*LFFT:
				break
			# Inner loop that actually reads the frames into the data array
			for j in range(framesWork):
				# Read in the next frame and anticipate any problems that could occur
				try:
					cFrame = tbw.readFrame(fh)
				except errors.eofError:
					break
				except errors.syncError:
					print "WARNING: Mark 5C sync error on frame #%i" % (int(fh.tell())/tbw.FrameSize-1)
					continue
				if not cFrame.header.isTBW():
					continue
				
				stand = cFrame.header.parseID()
				# In the current configuration, stands start at 1 and go up to 10.  So, we
				# can use this little trick to populate the data array
				aStand = 2*(stand-1)
				if cFrame.header.frameCount % 10000 == 0 and config['verbose']:
					print "%3i -> %3i  %6.3f  %5i  %i" % (stand, aStand, cFrame.getTime(), cFrame.header.frameCount, cFrame.data.timeTag)

				# Actually load the data.  x pol goes into the even numbers, y pol into the 
				# odd numbers
				count = cFrame.header.frameCount - 1
				data[aStand,   count*nSamples:(count+1)*nSamples] = cFrame.data.xy[0,:]
				data[aStand+1, count*nSamples:(count+1)*nSamples] = cFrame.data.xy[1,:]

			# Calculate the spectra for this block of data and then weight the results by 
			# the total number of frames read.  This is needed to keep the averages correct.
			# NB:  The weighting is the same for the x and y polarizations because of how 
			# the data are packed in TBW
			freq, tempSpec = fxc.SpecMaster(data, LFFT=LFFT, window=config['window'], verbose=config['verbose'])
			for stand in xrange(masterSpectra.shape[1]):
				masterSpectra[i,stand,:] = tempSpec[stand,:]

			# Compute the 1 ms average power and the data range within each 1 ms window
			subSize = 1960
			nSegments = data.shape[1] / subSize
			
			print "Computing average power and data range in %i-sample intervals" % subSize
			pb = ProgressBar(max=data.shape[0])
			avgPower = numpy.zeros((antpols, nSegments), dtype=numpy.float32)
			dataRange = numpy.zeros((antpols, nSegments, 3), dtype=numpy.int16)
			for s in xrange(data.shape[0]):
				for p in xrange(nSegments):
					subData = data[s,(p*subSize):((p+1)*subSize)]
					avgPower[s,p] = numpy.mean( numpy.abs(subData) )
					dataRange[s,p,0] = subData.min()
					dataRange[s,p,1] = subData.mean()
					dataRange[s,p,2] = subData.max()

					### This little block here looks for likely saturation events and save
					### the raw time series around them into individual NPZ files for stand
					### number 14.
					#if (dataRange[s,p,0] < -1000 or dataRange[s,p,0] > 1000) and antennas[s].stand.id == 14:
						#subData = data[s,((p-1)*1960):((p+2)*1960)]
						#satFileName = 'stand-14-pol-%i-%i.npz' % (antennas[s].pol, (p-1)*1960)
						#print satFileName
						#numpy.savez(satFileName, start=(p-1)*1960, data=subData)
				pb.inc(amount=1)
				if pb.amount != 0 and pb.amount % 10 == 0:
					sys.stdout.write(pb.show()+'\r')
					sys.stdout.flush()
			sys.stdout.write(pb.show()+'\r')
			sys.stdout.write('\n')
			sys.stdout.flush()

			# We don't really need the data array anymore, so delete it
			del(data)

		# Apply the cable loss corrections, if requested
		if config['applyGain']:
			for s in xrange(masterSpectra.shape[1]):
				currGain = antennas[s].cable.gain(freq)
				for c in xrange(masterSpectra.shape[0]):
					masterSpectra[c,s,:] /= currGain
					
		
		# Now that we have read through all of the chunks, perform the final averaging by
		# dividing by all of the chunks
		spec = masterSpectra.mean(axis=0)
		
		# Estimate the dipole resonance frequencies
		print "Computing dipole resonance frequencies"
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
			if pb.amount != 0 and pb.amount % 10 == 0:
				sys.stdout.write(pb.show()+'\r')
				sys.stdout.flush()
		sys.stdout.write(pb.show()+'\r')
		sys.stdout.write('\n')
		sys.stdout.flush()
		
		numpy.savez("%s.npz" % base, date=str(beginDate), freq=freq, masterSpectra=masterSpectra, resFreq=resFreq, 
					avgPower=avgPower, dataRange=dataRange, ssmifContents=ssmifContents)
	else:
		dataDict = numpy.load("%s.npz" % base)
		freq = dataDict['freq']
		masterSpectra = dataDict['masterSpectra']
		resFreq = dataDict['resFreq']

		# Now that we have read through all of the chunks, perform the final averaging by
		# dividing by all of the chunks
		spec = masterSpectra.mean(axis=0)
	
	# Create a good template spectra
	specTemplate = numpy.median(spec, axis=0)
	specDiff = numpy.zeros(spec.shape[0])
	toCompare = numpy.where( (freq>32e6) & (freq<50e6) )[0]
	print len(toCompare)
	for i in xrange(spec.shape[0]):
		specDiff[i] = (spec[i,toCompare] / specTemplate[toCompare]).mean()
	specDiff = numpy.where( specDiff < 2, specDiff, 2)
	
	# Get the station
	standPos = numpy.array([[ant.stand.x, ant.stand.y, ant.stand.z] for ant in antennas if ant.pol == 0])
	
	# Plots
	if config['verbose']:
		fig = plt.figure()
		ax1 = fig.add_subplot(1, 2, 1)
		ax1.scatter(standPos[:,0], standPos[:,1], c=specDiff[0::2], s=40.0, alpha=0.50)
		## Add the fence as a dashed line
		ax1.plot([-44.315, 72.150, 44.077, -72.543, -44.315], 
				[-72.522, -44.277, 72.191, 43.972, -72.522], linestyle='--', color='k')
		### Add the shelter
		#ax1.plot([55.863, 58.144, 58.062, 55.791, 55.863], 
				#[45.946, 45.999, 51.849, 51.838, 45.946], linestyle='-', color='k')
		## Set the limits to just zoom in on the main stations
		ax1.set_xlim([-75, 75])
		ax1.set_ylim([-75, 75])		
		
		ax2 = fig.add_subplot(1, 2, 2)
		ax2.plot(freq/1e6, numpy.log10(specTemplate)*10, alpha=0.50)
		
		print "RBW: %.1f Hz" % (freq[1]-freq[0])
		plt.show()


if __name__ == "__main__":
	main(sys.argv[1:])

