#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import ephem
import numpy

from lsl.common.stations import lwa1
from lsl.reader import tbn
from lsl.reader import errors
from lsl.reader.buffer import TBNFrameBuffer
from lsl.astro import unix_to_utcjd, DJD_OFFSET

from matplotlib import pyplot as plt


# List of bright radio sources and pulsars in PyEphem format
_srcs = ["ForA,f|J,03:22:41.70,-37:12:30.0,1",
         "TauA,f|J,05:34:32.00,+22:00:52.0,1", 
         "VirA,f|J,12:30:49.40,+12:23:28.0,1",
         "HerA,f|J,16:51:08.15,+04:59:33.3,1", 
         "SgrA,f|J,17:45:40.00,-29:00:28.0,1", 
         "CygA,f|J,19:59:28.30,+40:44:02.0,1", 
         "CasA,f|J,23:23:27.94,+58:48:42.4,1",]


def main(args):
	# The task at hand
	filename = args[0]
	
	# The station
	observer = lwa1.getObserver()
	antennas = lwa1.getAntennas()
	
	# The file's parameters
	fh = open(filename, 'rb')
	nFramesFile = os.path.getsize(filename) / tbn.FrameSize
	srate = tbn.getSampleRate(fh)
	antpols = len(antennas)
	
	# Reference antenna
	ref = 258
	for i,a in enumerate(antennas):
		if a.stand.id == ref and a.pol == 0:
			refX = i
		elif a.stand.id == ref and a.pol == 1:
			refY = i
		else:
			pass
	
	# Integration time (seconds and frames)
	tInt = 10.0
	nFrames = int(round(tInt*srate/512*antpols))
	tInt = nFrames / antpols * 512 / srate
	
	# Total run length
	nChunks = 1
	
	# Read in the first frame and get the date/time of the first sample 
	# of the frame.  This is needed to get the list of stands.
	junkFrame = tbn.readFrame(fh)
	fh.seek(-tbn.FrameSize, 1)
	startFC = junkFrame.header.frameCount
	centralFreq = junkFrame.getCentralFreq()
	beginDate = ephem.Date(unix_to_utcjd(junkFrame.getTime()) - DJD_OFFSET)
	
	observer.date = beginDate
	srcs = [ephem.Sun(),]
	for line in _srcs:
		srcs.append( ephem.readdb(line) )
	
	for i in xrange(len(srcs)):
		srcs[i].compute(observer)
		
		if srcs[i].alt > 0:
			print "source %s: alt %.1f degrees, az %.1f degrees" % (srcs[i].name, srcs[i].alt*180/numpy.pi, srcs[i].az*180/numpy.pi)

	# File summary
	print "Filename: %s" % filename
	print "Date of First Frame: %s" % str(beginDate)
	print "Ant/Pols: %i" % antpols
	print "Sample Rate: %i Hz" % srate
	print "Tuning Frequency: %.3f Hz" % centralFreq
	print "Frames: %i (%.3f s)" % (nFramesFile, 1.0 * nFramesFile / antpols * 512 / srate)
	print "---"
	print "Integration: %.3f s (%i frames; %i frames per stand/pol)" % (tInt, nFrames, nFrames / antpols)
	print "Chunks: %i" % nChunks
	
	junkFrame = tbn.readFrame(fh)
	while junkFrame.header.frameCount < startFC+3:
		junkFrame = tbn.readFrame(fh)
	fh.seek(-tbn.FrameSize, 1)
	
	# Create the FrameBuffer instance
	buffer = TBNFrameBuffer(stands=range(1,antpols/2+1), pols=[0, 1])
	
	# Create the phase average and times
	LFFT = 512
	times = numpy.zeros(nChunks, dtype=numpy.float64)
	fullVis   = numpy.zeros((nChunks, antpols, LFFT), dtype=numpy.complex64)
	simpleVis = numpy.zeros((nChunks, antpols), dtype=numpy.complex64)
	
	# Go!
	k = 0
	for i in xrange(nChunks):
		# Find out how many frames remain in the file.  If this number is larger
		# than the maximum of frames we can work with at a time (maxFrames),
		# only deal with that chunk
		framesRemaining = nFramesFile - k
		if framesRemaining > nFrames:
			framesWork = nFrames
			data = numpy.zeros((antpols, framesWork/antpols*512), dtype=numpy.complex64)
		else:
			framesWork = framesRemaining + antpols*buffer.nSegments
			data = numpy.zeros((antpols, framesWork/antpols*512), dtype=numpy.complex64)
		print "Working on chunk %i, %i frames remaining" % (i+1, framesRemaining)
		
		count = [0 for a in xrange(len(antennas))]
		
		j = 0
		fillsWork = framesWork / antpols
		# Inner loop that actually reads the frames into the data array
		while j < fillsWork:
			try:
				cFrame = tbn.readFrame(fh)
				k = k + 1
			except errors.eofError:
				break
			except errors.syncError:
				#print "WARNING: Mark 5C sync error on frame #%i" % (int(fh.tell())/tbn.FrameSize-1)
				continue
					
			buffer.append(cFrame)
			cFrames = buffer.get()

			if cFrames is None:
				continue
				
			for cFrame in cFrames:
				stand,pol = cFrame.header.parseID()
				
				# In the current configuration, stands start at 1 and go up to 260.  So, we
				# can use this little trick to populate the data array
				aStand = 2*(stand-1)+pol
				
				# Save the time
				if j == 0 and aStand == 0:
					times[i] = cFrame.getTime()
				
				data[aStand, count[aStand]*512:(count[aStand]+1)*512] = cFrame.data.iq
				
				# Update the counters so that we can average properly later on
				count[aStand] = count[aStand] + 1
			
			j += 1
			
		# Mask
		#bad = numpy.where( numpy.abs(data) >= 90 )
		#data[bad] *= 0.0
		
		# Simple correlation
		for l in xrange(520):
			if l % 2 == 0:
				simpleVis[i,l] = (data[l,:]*data[refX,:].conj()).mean()
			else:
				simpleVis[i,l] = (data[l,:]*data[refY,:].conj()).mean()
	
	# Save the data
	try:
		outname = filename.replace('.dat', '-vis.npz')
		numpy.savez(outname, ref=ref, refX=refX, refY=refY, tInt=tInt, centralFreq=centralFreq, times=times, fullVis=fullVis, simpleVis=simpleVis)
	except:
		pass
	
	# Plot the first 20
	fig = plt.figure()
	for i in xrange(20):
		ax = fig.add_subplot(5, 4, i+1)
		ax.plot(phase[:,i])
	plt.show()


if __name__ == "__main__":
	main(sys.argv[1:])
	
