#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
For a given MJD value or list of MJD values, return the range of local times 
associated with that MJD.

$Rev$
$LastChangedBy$
$LastChangedDate$
"""

import sys
import pytz
import getopt

from lsl.common.mcs import mjdmpm2datetime
from datetime import datetime

MST = pytz.timezone('US/Mountain')
UTC = pytz.utc


def usage(exitCode=None):
	print """mjd2local.py - For a given MJD value or list of MJD values, return
the range of local times associated with that MJD.

Usage: mjd2local.py [OPTIONS] MJD [MJD [MJD [...]]]

Options:
-h, --help                  Display this help information
-p, --pairs                 Interpret the input as MJD, MPM pairs
"""

	if exitCode is not None:
		sys.exit(exitCode)
	else:
		return True


def parseOptions(args):
	config = {}
	config['pairs'] = False

	# Read in and process the command line flags
	try:
		opts, args = getopt.getopt(args, "hp", ["help", "pairs"])
	except getopt.GetoptError, err:
		# Print help information and exit:
		print str(err) # will print something like "option -a not recognized"
		usage(exitCode=2)
	
	# Work through opts
	for opt, value in opts:
		if opt in ('-h', '--help'):
			usage(exitCode=0)
		if opt in ('-p', '--pairs'):
			config['pairs'] = True
		else:
			assert False
	
	# Add in arguments
	config['args'] = args

	# Return configuration
	return config


def main(args):
	config = parseOptions(args)

	if not config['pairs']:
		for arg in config['args']:
			mjd1 = int(arg)
			mjd2 = float(mjd1) + 0.99999

			d1 = mjdmpm2datetime(mjd1, 0)
			d1  = d1.astimezone(MST)

			d2 = mjdmpm2datetime(mjd2, 0)
			d2  = d2.astimezone(MST)

			print "MJD: %i" % mjd1
			print "Localtime: %s to %s" % (d1.strftime("%B %d, %Y at %H:%M:%S %Z"), d2.strftime("%B %d, %Y at %H:%M:%S %Z"))
	else:
		for arg in zip(config['args'][0::2], config['args'][1::2]):
			mjd, mpm = [int(i) for i in arg]
			d = mjdmpm2datetime(mjd, mpm)
			d = d.astimezone(MST)

			print "MJD: %i, MPM: %i" % (mjd, mpm)
			print "Localtime: %s" % d.strftime("%B %d, %Y at %H:%M:%S %Z")


if __name__ == "__main__":
	main(sys.argv[1:])




