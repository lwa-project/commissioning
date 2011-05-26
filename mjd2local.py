#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import pytz
import getopt

from lsl import astro
from datetime import datetime

MST = pytz.timezone('US/Mountain')
UTC = pytz.utc


def usage(exitCode=None):
	print """mjd2local.py - For a given MJD value or list of MJD values, return
the range of local times associated with that MJD.

Usage: mjd2local.py [OPTIONS] MJD [MJD [MJD [...]]]

Options:
-h, --help                  Display this help information
"""

	if exitCode is not None:
		sys.exit(exitCode)
	else:
		return True


def parseOptions(args):
	config = {}

	# Read in and process the command line flags
	try:
		opts, args = getopt.getopt(args, "h", ["help",])
	except getopt.GetoptError, err:
		# Print help information and exit:
		print str(err) # will print something like "option -a not recognized"
		usage(exitCode=2)
	
	# Work through opts
	for opt, value in opts:
		if opt in ('-h', '--help'):
			usage(exitCode=0)
		else:
			assert False
	
	# Add in arguments
	config['args'] = args

	# Return configuration
	return config


def main(args):
	config = parseOptions(args)

	for arg in config['args']:
		mjd1 = int(arg)
		mjd2 = float(mjd1) + 0.99999

		j1 = astro.mjd_to_jd(mjd1)
		t1  = astro.get_timet_from_julian(j1)
		d1  = UTC.localize(datetime.utcfromtimestamp(t1))
		d1  = d1.astimezone(MST)

		j2 = astro.mjd_to_jd(mjd2)
		t2  = astro.get_timet_from_julian(j2)
		d2  = UTC.localize(datetime.utcfromtimestamp(t2))
		d2  = d2.astimezone(MST)

		print "MJD: %i" % mjd1
		print "Localtime: %s to %s" % (d1.strftime("%B %d, %Y at %H:%M:%S %Z"), d2.strftime("%B %d, %Y at %H:%M:%S %Z"))


if __name__ == "__main__":
	main(sys.argv[1:])




