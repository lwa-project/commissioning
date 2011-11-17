#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Read in SSMIF file and create a set of DRX gain and delay files for a given 
frequency and topogentric pointing center.

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
	print """generateDelays.py - Read in a SSMIF file and create a set of DRX
gain and delay files for a given frequency and topogentric pointing center.

Usage: generateDelays.py [OPTIONS] file

Options:
-h, --help                  Display this help information
-f, --frequency             Frequency in MHz to calculate the gain/delays for 
                            (Default = 65 MHz)
-a, --azimuth               Azimuth east of north in degrees for the pointing center
                            (Default = 90 degrees)
-e, --elevation             Elevation above the horizon in degrees for the pointing
                            center (Default = 90 degrees)
-g, --gain                  DRX antenna gain (Default = 1.0000)
"""
	if exitCode is not None:
		sys.exit(exitCode)
	else:
		return True


def parseOptions(args):
	config = {}
	# Command line flags - default values
	config['freq'] = 65.0e6
	config['az'] = 90.0
	config['el'] = 90.0
	config['gain'] = 1.0000
	config['args'] = []

	# Read in and process the command line flags
	try:
		opts, args = getopt.getopt(args, "hf:a:e:g:", ["help", "frequency=", "azimuth=", "elevation=", "gain="])
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

	dftBase = 'beams_%iMHz_%iaz_%iel_%03ibg' % (config['freq']/1e6, config['az'], config['el'], config['gain']*100)
	gftBase = 'beams_%iMHz_%iaz_%iel_%03ibg' % (config['freq']/1e6, config['az'], config['el'], config['gain']*100)

	print "Calculating delays for az. %.2f, el. %.2f at %.2f MHz" % (config['az'], config['el'], config['freq']/1e6)
	delays = beamformer.calcDelay(antennas, freq=config['freq'], azimuth=config['az'], elevation=config['el'])
	delays *= 1e9
	delays = delays.max() - delays
	junk = delay.list2delayfile('.', dftBase, delays)

	print "Setting gains for %i good inputs, %i bad inputs" % (len(antennas)-len(bad), len(bad))
	bgain = config['gain']
        bgain_cross = 0.0000
        gains = [[bgain, bgain_cross, bgain_cross, bgain]]*260 # initialize gain list
	for d in digs[bad]:
		# Digitizers start at 1, list indicies at 0
		i = d - 1
		gains[i/2] = [0,0,0,0]
	junk = gain.list2gainfile('.', gftBase, gains)

	print "\nDelay and gain files are:\n %s.dft\n %s.gft" % (dftBase, gftBase)


if __name__ == "__main__":
	main(sys.argv[1:])

