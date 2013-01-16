#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Split multi-vis NPZ files that are generated from combined TBN observations into
single NPZ files, one for each frequency.

Usage:
./splitMultiVis.py <NPZ multi-vis. file> [<NPZ multi-vis. file> [...]]

$Rev$
$LastChangedBy$
$LastChangedDate$
"""

import os
import sys
import numpy
import getopt

from datetime import datetime


def usage(exitCode=None):
	print """splitMultiVis.py - Split multi-vis NPZ files that are generated from combined TBN
observations into single NPZ files, one for each frequency.

Usage: splitMultiVis.py [OPTIONS] file [file [...]]

Options:
-h, --help              Display this help information
-m, --min-integrations  Minimum number of integrations needed to keep
                        split out a frequency (default = 20)
"""

	if exitCode is not None:
		sys.exit(exitCode)
	else:
		return True


def parseOptions(args):
	config = {}
	# Command line flags - default values
	config['minInts'] = 20
	
	# Read in and process the command line flags
	try:
		opts, args = getopt.getopt(args, "hm:", ["help", "min-integrations="])
	except getopt.GetoptError, err:
		# Print help information and exit:
		print str(err) # will print something like "option -a not recognized"
		usage(exitCode=2)
	
	# Work through opts
	for opt, value in opts:
		if opt in ('-h', '--help'):
			usage(exitCode=0)
		elif opt in ('-m', '--metadata'):
			config['minInts'] = int(value)
		else:
			assert False
	
	# Add in arguments
	config['args'] = args

	# Return configuration
	return config


def main(args):
	config = parseOptions(args)
	filenames = config['args']

	for filename in filenames:
		dataDict = numpy.load(filename)
	
		print "Working on file '%s'" % filename

		# Load in the data
		refAnt = dataDict['ref'].item()
		refX   = dataDict['refX'].item()
		refY   = dataDict['refY'].item()
		tInt = dataDict['tInt'].item()
	
		times = dataDict['times']
		phase = dataDict['simpleVis']
	
		centralFreqs = dataDict['centralFreqs']

		# Find the unique sets of (non-zero) frequencies and report
		uFreq = numpy.unique(centralFreqs)
		uFreq = uFreq[numpy.where(uFreq != 0)]
		print "  Found %i unique frequencies from %.3f to %.3f MHz" % (len(uFreq), uFreq.min()/1e6, uFreq.max()/1e6)
	
		# Report on the start time
		beginDate = datetime.utcfromtimestamp(times[0])
		print "  Start date/time of data: %s UTC" % beginDate
	
		# Split
		for i,f in enumerate(uFreq):
			## Select what we need and trim off the last index to deal 
			## with frequency changes
			toKeep = numpy.where( centralFreqs == f )[0]
			toKeep = toKeep[:-1]
		
			## Sub-sets of `times` and `phase`
			subTimes = times[toKeep]
			subPhase = phase[toKeep,:]
		
			if len(toKeep) < config['minInts']:
				print "  -> Skipping %.3f MHz with only %i integrations (%.1f s)" % (f, len(toKeep), len(toKeep)*tInt)
				continue
		
			## Save the split data to its own file
			outname = os.path.splitext(filename)[0]
			outname = outname.replace('-multi', '')
			outname = "%s-%03i.npz" % (outname, i+1)
			print "  Saving visibility data for %.3f MHz to '%s'" % (f/1e6, outname)
			numpy.savez(outname, ref=refAnt, refX=refX, refY=refY, tInt=tInt, centralFreq=f, 
						times=subTimes, simpleVis=subPhase)


if __name__ == "__main__":
	main(sys.argv[1:])
	
