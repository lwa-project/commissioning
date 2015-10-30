#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Run through a DR spectrometer file and determine if it is bad or not.

$Rev$
$LastChangedBy$
$LastChangedDate$
"""

import os
import sys
import ephem
import numpy
import getopt
from datetime import datetime

from lsl import astro
from lsl.reader import drx, drspec, errors


def usage(exitCode=None):
	print """drspecFileCheck.py - Run through a DR spectrometer file and determine if it is bad or not.

Usage: drspecFileCheck.py [OPTIONS] filename

Options:
-h, --help         Display this help information
-l, --length       Length of time in seconds to analyze (default 1 s)
-s, --skip         Skip period in seconds between chunks (default 900 s)
-t, --trim-level   Trim level for power analysis with clipping (default 49)
"""

	if exitCode is not None:
		sys.exit(exitCode)
	else:
		return None


def parseOptions(args):
	config = {}
	config['length'] = 1.0
	config['skip'] = 900.0
	config['trim'] = 49
	
	try:
		opts, args = getopt.getopt(args, "hl:s:t:", ["help", "length=", "skip=", "trim-level="])
	except getopt.GetoptError, err:
		# Print help information and exit:
		print str(err) # will print something like "option -a not recognized"
		usage(exitCode=2)

	# Work through opts
	for opt, value in opts:
		if opt in ('-h', '--help'):
			usage(exitCode=0)
		elif opt in ('-l', '--length'):
			config['length'] = float(value)
		elif opt in ('-s', '--skip'):
			config['skip'] = float(value)
		elif opt in ('-t', '--trim-level'):
			config['trim'] = int(value)
		else:
			assert False
			
	# Add in arguments
	config['args'] = args
	
	# Return
	return config


def main(args):
	config = parseOptions(args)
	filename = config['args'][0]
	
	fh = open(filename, "rb")
	
	try:
		for i in xrange(5):
			junkFrame = drx.readFrame(fh)
		raise RuntimeError("ERROR: '%s' appears to be a raw DRX file, not a DR spectrometer file" % filename)
	except errors.syncError:
		fh.seek(0)
		
	# Interrogate the file to figure out what frames sizes to expect, now many 
	# frames there are, and what the transform length is
	FrameSize = drspec.getFrameSize(fh)
	nFrames = os.path.getsize(filename) / FrameSize
	nChunks = nFrames
	LFFT = drspec.getTransformSize(fh)
	
	# Read in the first frame to figure out the DP information
	junkFrame = drspec.readFrame(fh)
	fh.seek(-FrameSize, 1)
	srate = junkFrame.getSampleRate()
	t0 = junkFrame.getTime()
	tInt = junkFrame.header.nInts*LFFT/srate
	
	# Update the file contents
	beam = junkFrame.parseID()
	centralFreq1 = junkFrame.getCentralFreq(1)
	centralFreq2 = junkFrame.getCentralFreq(2)
	srate = junkFrame.getSampleRate()
	dataProducts = junkFrame.getDataProducts()
	t0 = junkFrame.getTime()
	tInt = junkFrame.header.nInts*LFFT/srate
	beginDate = datetime.utcfromtimestamp(junkFrame.getTime())
        
	# Report
	print "Filename: %s" % filename
	print "Date of First Frame: %s" % beginDate
	print "Beam: %i" % beam
	print "Sample Rate: %i Hz" % srate
	print "Tuning Frequency: %.3f Hz (1); %.3f Hz (2)" % (centralFreq1, centralFreq2)
	print "Data Products: %s" % ','.join(dataProducts)
	print "Frames: %i (%.3f s)" % (nFrames, nFrames*tInt)
	print "---"
	print "Transform Length: %i" % LFFT
	print "Integration: %.3f s" % tInt
	
	# Convert chunk length to total frame count
	chunkLength = int(config['length'] / tInt)
	
	# Convert chunk skip to total frame count
	chunkSkip = int(config['skip'] / tInt)
	
	# Output arrays
	clipFraction = []
	meanPower = []
	
	# Go!
	i = 1
	done = False
	print "   |%sClipping%s |%sPower %s |" % (" "*(8*len(dataProducts)-4), " "*(8*len(dataProducts)-4), " "*(6*len(dataProducts)-3), " "*(6*len(dataProducts)-3))
	out = "   |      1X      1Y      2X      2Y |"
	for t in (1, 2):
		for dp in dataProducts:
			out += "%6s" % ("%i%s" % (t, dp))
	out += " |"
	print out
	print "-"*len(out)
	
	while True:
		count = {0:0, 1:0, 2:0, 3:0}
		sats = numpy.empty((4,chunkLength), dtype=numpy.float32)
		data = numpy.empty((2*len(dataProducts),chunkLength*LFFT), dtype=numpy.float32)
		for j in xrange(chunkLength):
			# Read in the next frame and anticipate any problems that could occur
			try:
				cFrame = frame = drspec.readFrame(fh)
			except errors.eofError:
				done = True
				break
			except errors.syncError:
				continue
				
			for t in (1,2):
				for p,dp in enumerate(dataProducts):
					l = len(dataProducts)*(t-1) + p
					data[l,j*LFFT:(j+1)*LFFT] = getattr(cFrame.data, '%s%i' % (dp,t-1))
			sats[:,j] = numpy.array(cFrame.data.saturations) / (tInt*srate)
					
		if done:
			break
			
		else:
			clipFraction.append( sats.mean(axis=1) )
			meanPower.append( data.mean(axis=1) )
			
			clip = clipFraction[-1]
			power = meanPower[-1]
			
			out = "%2i | %6.2f%% %6.2f%% %6.2f%% %6.2f%% |" % (i, clip[0]*100.0, clip[1]*100.0, clip[2]*100.0, clip[3]*100.0)
			for t in (1, 2):
				for p in xrange(len(dataProducts)):
					out += " %5.2f" % (power[len(dataProducts)*(t-1)+p],)
			out += " |"
			print out
		
			i += 1
			fh.seek(FrameSize*chunkSkip, 1)
			
	clipFraction = numpy.array(clipFraction)
	meanPower = numpy.array(meanPower)
	
	clip = clipFraction.mean(axis=0)
	power = meanPower.mean(axis=0)
	
	print "-"*len(out)
	out = "%2s | %6.2f%% %6.2f%% %6.2f%% %6.2f%% |" % ('M', clip[0]*100.0, clip[1]*100.0, clip[2]*100.0, clip[3]*100.0)
	for t in (1, 2):
		for p in xrange(len(dataProducts)):
			out += " %5.2f" % (power[len(dataProducts)*(t-1)+p],)
	out += " |"
	print out


if __name__ == "__main__":
	main(sys.argv[1:])
	
