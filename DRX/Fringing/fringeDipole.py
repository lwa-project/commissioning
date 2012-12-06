#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Script to fringe special DRX files that have one dipole on X pol. and another
dipole on Y pol.  The visibilites are written to a NPZ file.

$Rev$
$LastChangedBy$
$LastChangedDate$
"""

import os
import sys
import numpy
import getopt

from lsl.reader import drx
from lsl.common.dp import fS
from lsl.common import stations
from lsl.astro import unix_to_utcjd, DJD_OFFSET
from lsl.correlator import fx as fxc
from lsl.common.progress import ProgressBar

from matplotlib import pyplot as plt


def usage(exitCode=None):
	print """fringeDipole.py - Take a DRX file with a dipole on X pol. and a dipole
on the Y pol. and cross correlate it.

Usage: fringeDipole.py [OPTION] <dipole_ID_X> <dipole_ID_Y> <DRX_file>

Options:
-h, --help                  Display this help information
-l, --fft-length            Set FFT length (default = 512)
-t, --avg-time              Window to average visibilities in time (seconds; 
                            default = 4 s)
"""

	if exitCode is not None:
		sys.exit(exitCode)
	else:
		return True


def parseOptions(args):
	config = {}
	config['avgTime'] = 4.0
        config['LFFT'] = 512

	# Read in and process the command line flags
	try:
		opts, args = getopt.getopt(args, "hl:t:", ["help", "fft-length=", "avg-time="])
	except getopt.GetoptError, err:
		# Print help information and exit:
		print str(err) # will print something like "option -a not recognized"
		usage(exitCode=2)
	
	# Work through opts
	for opt, value in opts:
		if opt in ('-h', '--help'):
			usage(exitCode=0)
		elif opt in ('-l', '--fft-length'):
                        config['LFFT'] = int(value)
                elif opt in ('-t', '--avg-time'):
                        config['avgTime'] = float(value)
		else:
			assert False
	
	# Add in arguments
	config['args'] = args

	# Return configuration
	return config


def main(args):
	config = parseOptions(args)

	LFFT = config['LFFT']

	stand1 = int(config['args'][0])
	stand2 = int(config['args'][1])
	filenames = config['args'][2:]
	
	# Build up the station
	site = stations.lwa1
	
	# Figure out which antennas we need
	antennas = []
	for ant in site.getAntennas():
		if ant.stand.id == stand1 and ant.pol == 0:
			antennas.append(ant)
	for ant in site.getAntennas():
		if ant.stand.id == stand2 and ant.pol == 0:
			antennas.append(ant)
	
	# Loop through the input files...
	for filename in filenames:
		fh = open(filename, "rb")
		nFramesFile = os.path.getsize(filename) / drx.FrameSize
		junkFrame = drx.readFrame(fh)
		fh.seek(0)
	
		beam, tune, pol = junkFrame.parseID()
		srate = junkFrame.getSampleRate()
	
		tunepols = drx.getFramesPerObs(fh)
		tunepols = tunepols[0] + tunepols[1] + tunepols[2] + tunepols[3]
		beampols = tunepols
		
		tnom = junkFrame.header.timeOffset
		tStart = junkFrame.getTime()
	
		# Get the DRX frequencies
		cFreq1 = 0.0
		cFreq2 = 0.0
		for i in xrange(4):
			junkFrame = drx.readFrame(fh)
			b,t,p = junkFrame.parseID()
			if p == 0 and t == 1:
				try:
					cFreq1 = junkFrame.getCentralFreq()
				except AttributeError:
					from lsl.common.dp import fS
					cFreq1 = fS * ((junkFrame.data.flags>>32) & (2**32-1)) / 2**32
			elif p == 0 and t == 2:
				try:
					cFreq2 = junkFrame.getCentralFreq()
				except AttributeError:
					from lsl.common.dp import fS
					cFreq2 = fS * ((junkFrame.data.flags>>32) & (2**32-1)) / 2**32
			else:
				pass
		fh.seek(-4*drx.FrameSize, 1)
	
		# Align the files as close as possible by the time tags and then make sure that
		# the first frame processed is from tuning 1, pol 0.
		junkFrame = drx.readFrame(fh)
		beam, tune, pol = junkFrame.parseID()
		pair = 2*(tune-1) + pol
		j = 0
		while pair != 0:
			junkFrame = drx.readFrame(fh)
			beam, tune, pol = junkFrame.parseID()
			pair = 2*(tune-1) + pol
			j += 1
		fh.seek(-drx.FrameSize, 1)
		print "Shifted beam %i data by %i frames (%.4f s)" % (beam, j, j*4096/srate/4)
	
		# Set integration time
		tInt = config['avgTime']
		nFrames = int(round(tInt*srate/4096))
		tInt = nFrames*4096/srate
	
		# Read in some data
		tFile = nFramesFile / 4 * 4096 / srate
	
		# Report
		print "Filename: %s" % filename
		print "  Sample Rate: %i Hz" % srate
		print "  Tuning 1: %.1f Hz" % cFreq1
		print "  Tuning 2: %.1f Hz" % cFreq2
		print "  ==="
		print "  Integration Time: %.3f s" % tInt
		print "  Integrations in File: %i" % int(tFile/tInt)

		nChunks = int(tFile/tInt)
		pb = ProgressBar(max=nChunks)
		for i in xrange(nChunks):
			junkFrame = drx.readFrame(fh)
			tStart = junkFrame.getTime()
			fh.seek(-drx.FrameSize, 1)

			count1 = [0,0]
			data1 = numpy.zeros((2, 4096*nFrames), dtype=numpy.complex64)
			count2 = [0,0]
			data2 = numpy.zeros((2, 4096*nFrames), dtype=numpy.complex64)
			for j in xrange(nFrames):
				for k in xrange(4):
					cFrame = drx.readFrame(fh)
					beam, tune, pol = cFrame.parseID()
					pair = 2*(tune-1) + pol

					if tune == 1:
						data1[pol, count1[pol]*4096:(count1[pol]+1)*4096] = cFrame.data.iq
						count1[pol] += 1
					else:
						data2[pol, count2[pol]*4096:(count2[pol]+1)*4096] = cFrame.data.iq
						count2[pol] += 1
					
			# Correlate
			blList1, freq1, vis1 = fxc.FXMaster(data1, antennas, LFFT=LFFT, Overlap=1, IncludeAuto=True, verbose=False, SampleRate=srate, CentralFreq=cFreq1, Pol='XX', ReturnBaselines=True, GainCorrect=False, ClipLevel=0)
		
			blList2, freq2, vis2 = fxc.FXMaster(data2, antennas, LFFT=LFFT, Overlap=1, IncludeAuto=True, verbose=False, SampleRate=srate, CentralFreq=cFreq2, Pol='XX', ReturnBaselines=True, GainCorrect=False, ClipLevel=0)
	
			if nChunks != 1:
				outfile = os.path.splitext(filename)[0]
				outfile = "%s-vis-%04i.npz" % (outfile, i+1)
			else:
				outfile = os.path.splitext(filename)[0]
				outfile = "%s-vis.npz" % outfile
			numpy.savez(outfile, srate=srate, freq1=freq1, vis1=vis1, freq2=freq2, vis2=vis2, tStart=tStart, tInt=tInt, stands=numpy.array([stand1, stand2]))

			del data1
			del data2

			pb.inc(amount=1)
			sys.stdout.write(pb.show()+'\r')
			sys.stdout.flush()

		sys.stdout.write(pb.show()+'\r')
		sys.stdout.write('\n')
		sys.stdout.flush()


		# Plot
		fig = plt.figure()
		i = 0
		for bl, vi in zip(blList1, vis1):
			ax = fig.add_subplot(4, 3, i+1)
			ax.plot(freq1/1e6, numpy.unwrap(numpy.angle(vi)))
			ax.set_title('Stand %i - Stand %i' % (bl[0].stand.id, bl[1].stand.id))
			ax = fig.add_subplot(4, 3, i+4)
			ax.plot(freq1/1e6, numpy.abs(vi))
			i += 1

			coeff = numpy.polyfit(freq1, numpy.unwrap(numpy.angle(vi)), 1)
			#print coeff[0]/2/numpy.pi*1e9, coeff[1]*180/numpy.pi
		
		i = 6
		for bl, vi in zip(blList2, vis2):
			ax = fig.add_subplot(4, 3, i+1)
			ax.plot(freq2/1e6, numpy.unwrap(numpy.angle(vi)))
			ax.set_title('Stand %i - Stand %i' % (bl[0].stand.id, bl[1].stand.id))
			ax = fig.add_subplot(4, 3, i+4)
			ax.plot(freq2/1e6, numpy.abs(vi))
			i += 1

			coeff = numpy.polyfit(freq2, numpy.unwrap(numpy.angle(vi)), 1)
			#print coeff[0]/2/numpy.pi*1e9, coeff[1]*180/numpy.pi

		#plt.show()
	

if __name__ == "__main__":
	main(sys.argv[1:])
	
