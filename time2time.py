#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import math
import pytz
import datetime
import getopt


def usage(exitCode=None):
	print """time2time.py - Convert a local date/time string in the format of 
"YYYY-MM-DD HH:MM:SS[.SSS]" into MJD and MPM UTC values.  If no date/time string
is supplied, the current local date/time is used.

Usage: time2time.py [OPTIONS] YYYY-MM-DD HH:MM:SS[.SSS]

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


def time2time(dt=datetime.datetime.utcnow()):
	"""
	Convert a UTC datetime object into MJD and MPM.  Based on the JPL get_time()
	function found in DRX/test_analog_in2.py.
	"""
	
	year        = dt.year             
	month       = dt.month      
	day         = dt.day    

	hour        = dt.hour
	minute      = dt.minute
	second      = dt.second     
	millisecond = dt.microsecond / 1000

	# compute MJD         
	# adapted from http://paste.lisp.org/display/73536
	# can check result using http://www.csgnetwork.com/julianmodifdateconv.html
	a = (14 - month) // 12
	y = year + 4800 - a          
	m = month + (12 * a) - 3                    
	p = day + (((153 * m) + 2) // 5) + (365 * y)   
	q = (y // 4) - (y // 100) + (y // 400) - 32045
	mjdi = int(math.floor( (p+q) - 2400000.5))
	mjd = mjdi

	# compute MPM
	mpmi = int(math.floor( (hour*3600 + minute*60 + second)*1000 + millisecond ))
	mpm = mpmi
	return (mjd, mpm)


def main(args):
	config = parseOptions(args)
	
	MST = pytz.timezone('US/Mountain')
	UTC = pytz.utc
	
	if len(config['args']) == 0:
		dt = datetime.datetime.utcnow()
		dt = UTC.localize(dt)
	else:
		year, month, day = args[0].split('-', 2)
		hour, minute, second = args[1].split(':', 2)
		iSeconds = int(float(second))
		mSeconds = int(round((float(second) - iSeconds)*1000000))

		dt = MST.localize(datetime.datetime(int(year), int(month), int(day), int(hour), int(minute), iSeconds, mSeconds))
		dt = dt.astimezone(UTC)
	
	mjd, mpm = time2time(dt)
	
	print "Localtime: %s" % dt.astimezone(MST).strftime("%B %d, %Y at %H:%M:%S %Z")
	print "MJD: %i" % mjd
	print "MPM: %i" % mpm
	

if __name__ == "__main__":
	main(sys.argv[1:])

