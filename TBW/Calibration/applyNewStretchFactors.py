#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
$Rev$
$LastChangedBy$
$LastChangedDate$
"""

import os
import sys
import numpy
import getopt

from lsl.common.stations import parseSSMIF

# List of stands *not* to update
EXCLUDED_STANDS = []


def usage(exitCode=None):
	print """applyNewStretchFactors.py - Given an existing SSMIF and new stretch factors, build 
a new SSMIF.

Usage: applyNewStretchFactors.py [OPTIONS] SSMIF stretchFile

Options:
-h, --help             Display this help information
-e, --exclude          Comma seperated list of stands not to update
                       (default = update all stands)
-o, --output           Write output to the specified filename 
                       (default = write to screen)
"""
	
	if exitCode is not None:
		sys.exit(exitCode)
	else:
		return True


def parseConfig(args):
	config = {}
	# Command line flags - default values
	config['exclude'] = []
	config['outname'] = None
	config['showPlots'] = False

	# Read in and process the command line flags
	try:
		opts, arg = getopt.getopt(args, "he:o:", ["help", "exclude=", "output="])
	except getopt.GetoptError, err:
		# Print help information and exit:
		print str(err) # will print something like "option -a not recognized"
		usage(exitCode=2)
		
	# Work through opts
	for opt, value in opts:
		if opt in ('-h', '--help'):
			usage(exitCode=0)
		elif opt in ('-e', '--exclude'):
			config['exclude'] = [int(v, 10) for v in value.split(',')]
		elif opt in ('-o', '--output'):
			config['outname'] = value
		else:
			assert False
			
	# Add in arguments
	config['args'] = arg
	
	# Validate the input
	if len(config['args']) != 2:
		raise RuntimeError("Must provide both a SSMIF and stretch file")
		
	# Return configuration
	return config


def main(args):
	#
	# Parse the command line
	#
	config = parseConfig(args)
	
	#
	# Load in the data
	#
	ssmifContents = open(config['args'][0], 'r').readlines()
	site     = parseSSMIF(config['args'][0])
	dataFile = numpy.loadtxt(config['args'][1])
	
	#
	# Gather the station meta-data from its various sources
	#
	observer = site.getObserver()
	antennas = site.getAntennas()
	
	#
	# Match the new stretch factors to the antennas
	#
	factors = [1.0 for i in xrange(len(antennas))]
	for i in xrange(dataFile.shape[0]):
		dig, stretch, addDelay, rms, chi2 = dataFile[i,:]
		dig = int(dig)
		antenna = antennas[dig-1]
		if antenna.stand.id in config['exclude']:
			continue
			
		factors[antenna.id-1] = stretch
		
	#
	# Final results
	#
	if config['outname'] is not None:
		fh = open(config['outname'], 'w')
	else:
		fh = sys.stdout
		
	for line in ssmifContents:
		if line[0:8] == 'RPD_STR[':
			start = line.find('[')
			stop  = line.find(']')
			try:
				junk, toSave = line.split('#', 1)
				toSave = " # %s" % toSave
			except ValueError:
				toSave = "\n"
			
			antID = int(line[start+1:stop])
			fh.write("RPD_STR[%i]  %.4f%s" % (antID, factors[antID-1], toSave))
		else:
			fh.write(line)
			
	if config['outname'] is not None:
		fh.close()


if __name__ == "__main__":
	main(sys.argv[1:])
	