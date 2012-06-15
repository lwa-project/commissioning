#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Convert a local date/time string in the format of "YYYY/MM/DD HH:MM:SS[.SSS]" into 
MJD and MPM UTC values.  If no date/time string is supplied, the current local 
date/time is used.

$Rev$
$LastChangedBy$
$LastChangedDate$
"""

import sys
import math
import pytz
import datetime
import getopt

from lsl.common.mcs import datetime2mjdmpm


def usage(exitCode=None):
	print """time2time.py - Convert a local date/time string in the format of 
"YYYY-MM-DD HH:MM:SS[.SSS]" into MJD and MPM UTC values.  If no date/time string
is supplied, the current local date/time is used.

Usage: time2time.py [OPTIONS] YYYY/MM/DD HH:MM:SS[.SSS]

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
	
	MST = pytz.timezone('US/Mountain')
	UTC = pytz.utc
	
	if len(config['args']) == 0:
		dt = datetime.datetime.utcnow()
		dt = UTC.localize(dt)
	else:
		args[0] = args[0].replace('-', '/')
		year, month, day = args[0].split('/', 2)
		hour, minute, second = args[1].split(':', 2)
		iSeconds = int(float(second))
		mSeconds = int(round((float(second) - iSeconds)*1000000))

		dt = MST.localize(datetime.datetime(int(year), int(month), int(day), int(hour), int(minute), iSeconds, mSeconds))
		dt = dt.astimezone(UTC)
	
	mjd, mpm = datetime2mjdmpm(dt)
	
	print "Localtime: %s" % dt.astimezone(MST).strftime("%B %d, %Y at %H:%M:%S %Z")
	print "MJD: %i" % mjd
	print "MPM: %i" % mpm
	

if __name__ == "__main__":
	main(sys.argv[1:])

