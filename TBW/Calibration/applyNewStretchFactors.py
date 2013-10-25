#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import numpy
import tempfile

from lsl.common.stations import parseSSMIF, lwa1

# List of stands *not* to update
EXCLUDED_STANDS = []


def main(args):
	#
	# Load in the data
	#
	ssmifContents = open(args[0], 'r').readlines()
	site     = parseSSMIF(args[0])
	dataFile = numpy.loadtxt(args[1])
	
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
		if antenna.stand.id in EXCLUDED_STANDS:
			continue
			
		factors[antenna.id-1] = stretch
		
	#
	# Final results
	#
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
			print "RPD_STR[%i]  %.4f%s" % (antID, factors[antID-1], toSave),
		else:
			print line,


if __name__ == "__main__":
	main(sys.argv[1:])
	
