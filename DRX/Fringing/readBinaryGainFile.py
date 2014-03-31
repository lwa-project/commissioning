#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Simple script to read in a MCS binary packed DP delay file (.df) and print 
out the delays in ns.

$Rev$
$LastChangedBy$
$LastChangedDate$
"""

import os
import sys
import numpy
import struct

from lsl.common.dp import DPGtogain
from lsl.common.stations import lwa1


def main(args):
	filename = args[0]

	# Read in the entire file
	fh = open(filename, 'rb')
	data = fh.read()
	fh.close()
	
	# Unpack the raw delays (260*4 signed short ints)
	rawGains = struct.unpack('<1040h', data)

	# Convert to delays in ns
	gains = [DPGtogain(g) for g in rawGains]
	
	# Report
	ants = lwa1.getAntennas()[0::2]
	
	print "Std   X->x   X->y    Y->x   Y->y"
	print "----------------------------------"
	for i in xrange(len(ants)):
		xofx, xofy, yofx, yofy = gains[4*i:4*(i+1)]
		print "%3i   %6.4f %6.4f  %6.4f %6.4f" % (ants[i].stand.id, xofx, xofy, yofx, yofy)


if __name__ == "__main__":
	main(sys.argv[1:])

