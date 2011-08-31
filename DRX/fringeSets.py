#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Read in SSMIF file and create a set of DRX gain files (zero-delay for dipoles-only) 
that puts a single dipole or beam on the X pol and the outlier on the other.  The 
gain files are constructed such that all data is from X pol.

$Rev$
$LastChangedBy$
$LastChangedDate$
"""

import sys
import numpy
import getopt

import gain
import delay
from lsl.misc import beamformer
from lsl.common.stations import parseSSMIF

def usage(exitCode=None):
	print """fringeSets.py - Read in SSMIF file and create a set of DRX
gain files (zero-delay) that puts a single dipole or beam on the X pol and
the outlier on the other.

Usage: fringeSets.py [OPTIONS] file

Options:
-h, --help                  Display this help information
-f, --frequency             Beam only, frequency in MHz to calculate the gain/delays 
                            for (Default = 74 MHz)
-a, --azimuth               Beam only, azimuth east of north in degrees for the 
                            pointing center (Default = 90 degrees)
-e, --elevation             Beam only, elevation above the horizon in degrees for 
                            the pointing center (Default = 90 degrees)
-d, --dipole                Using a dipole instead of the beam (Default = use beam)
-r, --reference             Reference for the fringing (Default = stand #258)
-g, --gain                  DRX antenna gain (Default = 1.0000)
"""
	if exitCode is not None:
		sys.exit(exitCode)
	else:
		return True


def parseOptions(args):
	config = {}
	# Command line flags - default values
	config['beam'] = True
	config['freq'] = 74.0e6
	config['az'] = 90.0
	config['el'] = 90.0
	config['gain'] = 1.0000
	config['dipole'] = 1
	config['ref'] = 258
	config['args'] = []
	
	# Read in and process the command line flags
	try:
		opts, args = getopt.getopt(args, "hf:a:e:g:d:r:", ["help", "frequency=", "azimuth=", "elevation=", "gain=", "dipole=", "reference="])
	except getopt.GetoptError, err:
		# Print help information and exit:
		print str(err) # will print something like "option -a not recognized"
		usage(exitCode=2)
	
	# Work through opts
	for opt, value in opts:
		if opt in ('-h', '--help'):
			usage(exitCode=0)
		elif opt in ('-f', '--frequency'):
			config['freq'] = float(value)*1e6
		elif opt in ('-a', '--azimuth'):
			config['az'] = float(value)
		elif opt in ('-e', '--elevation'):
			config['el'] = float(value)
		elif opt in ('-g', '--gain'):
			config['gain'] = float(value)
		elif opt in ('-d', '--dipole'):
			config['beam'] = False
			config['dipole'] = int(value)
		elif opt in ('-r', '--reference'):
			config['ref'] = int(value)
		else:
			assert False
	
	# Add in arguments
	config['args'] = args

	# Return configuration
	return config


def main(args):
	config = parseOptions(args)
	filename = config['args'][0]

	station = parseSSMIF(filename)
	antennas = station.getAntennas()

	digs    = numpy.array([ant.digitizer  for ant in antennas])
	ants    = numpy.array([ant.id         for ant in antennas])
	stands  = numpy.array([ant.stand.id   for ant in antennas])
	pols    = numpy.array([ant.pol        for ant in antennas])
	antStat = numpy.array([ant.status     for ant in antennas])
	feeStat = numpy.array([ant.fee.status for ant in antennas])

	badStands = numpy.where( antStat != 3 )[0]
	badFees   = numpy.where( feeStat != 3 )[0]
	bad = numpy.where( (stands > 256) | (antStat != 3) | (feeStat != 3) )[0]
	print "Number of bad stands:   %3i" % len(badStands)
	print "Number of bad FEEs:     %3i" % len(badFees)
	print "---------------------------"
	print "Total number bad inuts: %3i" % len(bad)
	print " "
	
	if config['beam']:
		dftBase = 'fringe_%iMHz_%iaz_%iel_%iref' % (config['freq']/1e6, config['az'], config['el'], config['ref'])
		gftBase = 'fringe_%iMHz_%iaz_%iel_%iref' % (config['freq']/1e6, config['az'], config['el'], config['ref'])
	else:
		dftBase = None
		gftBase = 'fringe_%istand_%iref' % (config['dipole'], config['ref'])
	
	bgain = config['gain']
	if config['beam']:
		print "Calculating delays for az. %.2f, el. %.2f at %.2f MHz" % (config['az'], config['el'], config['freq']/1e6)
		delays = beamformer.calcDelay(antennas, freq=config['freq'], azimuth=config['az'], elevation=config['el'])
		delays *= 1e9
		delays[bad] = 0
		junk = delay.list2delayfile('.', dftBase, delays)
		
		print "Setting gains for %i good inputs, %i bad inputs" % (len(antennas)-len(bad), len(bad))
		
		gains = [[bgain, 0.0000, 0.0000, 0.0000]]*260 # initialize gain list
		for d in digs[bad]:
			# Digitizers start at 1, list indicies at 0
			i = d - 1
			gains[i/2] = [0.0000,]*4
			
		for i in xrange(len(stands)/2):
			# Put the reference stand in there all by itself
			if stands[2*i] == config['ref']:
				gains[i] = [0.0000, bgain, 0.0000, 0.0000]
	else:
		print "Setting gains for dipoles %i and %i" % (config['dipole'], config['ref'])
		
		gains = [[0.0000, 0.0000, 0.0000, 0.0000]]*260 # initialize gain list
		for i in xrange(len(stands)/2):
			# Put the fringing stand in there all by itself
			if stands[2*i] == config['dipole']:
				gains[i] = [bgain, 0.0000, 0.0000, 0.0000]
			
			# Put the reference stand in there all by itself
			if stands[2*i] == config['ref']:
				gains[i] = [0.0000, bgain, 0.0000, 0.0000]

	#for g,s,p in zip(gains, stands[::2], pols[::2]):
		#print g, s, p
	junk = gain.list2gainfile('.', gftBase, gains)

	if config['beam']:
		print "\nDelay and gain files are:\n %s.dft\n %s.gft" % (dftBase, gftBase)
	else:
		print "Gain file is:\n %s.gft" % gftBase


if __name__ == "__main__":
	main(sys.argv[1:])
