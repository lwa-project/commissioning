#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
New take on astroevents using PyEphem for calculations.  It can also take a
date in the form of YYYY/MM/DD from the command line to use a as base for 
its calculations.

$Rev$
$LastChangedBy$
$LastChangedDate$
"""

import os
import sys
import pytz
import math
import ephem
import getopt
from datetime import datetime

from lsl.common.stations import lwa1


# Time zones
_UTC = pytz.utc
_MST = pytz.timezone('US/Mountain')


# List of bright radio sources and pulsars in PyEphem format
_srcs = ["ForA,f|J,03:22:41.70,-37:12:30.0,1",
         "TauA,f|J,05:34:32.00,+22:00:52.0,1", 
         "VirA,f|J,12:30:49.40,+12:23:28.0,1",
         "HerA,f|J,16:51:08.15,+04:59:33.3,1", 
         "SgrA,f|J,17:45:40.00,-29:00:28.0,1", 
         "CygA,f|J,19:59:28.30,+40:44:02.0,1", 
         "CasA,f|J,23:23:27.94,+58:48:42.4,1",
	 "B0329+54,f|J,03:32:59.37,+54:34:43.6,1",
	 "B0809+74,f|J,08:14:59.44,+74:29:05.8,1", 
         "B0950+08,f|J,09:53:09.31,+07:55:35.8,1",
         "B1133+16,f|J,11:36:03.25,+15:51:04.5,1",
         "B1919+21,f|J,19:21:44.80,+21:53:01.8,1",
	 "J2145-0750,f|J,21:45:50.47,-07:50:18.3,1",]

def usage(exitCode=None):
	print """astroevents2.py - New take on the astroevents.py script included in LSL 0.5.0+
that uses PyEphem for its calculations and can make calculations for a different day.

Usage: astroevents2.py [OPTIONS] [YYYY/MM/DD [HH:MM:SS]]

Options:
-h, --help                  Display this help information
-p, --position-mode         Display the azimuth and elevation of sources above the
                            horizon.
"""

	if exitCode is not None:
		sys.exit(exitCode)
	else:
		return True


def parseOptions(args):
	config = {}
	config['positionMode'] = False

	# Read in and process the command line flags
	try:
		opts, args = getopt.getopt(args, "hp", ["help", "position-mode"])
	except getopt.GetoptError, err:
		# Print help information and exit:
		print str(err) # will print something like "option -a not recognized"
		usage(exitCode=2)
	
	# Work through opts
	for opt, value in opts:
		if opt in ('-h', '--help'):
			usage(exitCode=0)
		elif opt in ('-p', '--position-mode'):
			config['positionMode'] = True
		else:
			assert False
	
	# Add in arguments
	config['args'] = args

	# Return configuration
	return config


def main(args):
	# Parse the command line
	config = parseOptions(args)

	# Get LWA-1
	observer = lwa1.getObserver()
	print "Current site is %s at lat %s, lon %s" % (lwa1.name, observer.lat, observer.long)
	
	# Set the current time so we can find the "next" transit.  Go ahead
	# and report the time and the current LST (for reference)
	if len(config['args']) == 2:
		year, month, day = config['args'][0].split('/', 2)
		year = int(year)
		month = int(month)
		day = int(day)
		
		hour, minute, second = config['args'][1].split(':', 2)
		hour = int(hour)
		minute = int(minute)
		second = int(second)
		
		tNow = _MST.localize(datetime(year, month, day, hour, minute, second))
		tNow = tNow.astimezone(_UTC)
	
	elif len(config['args']) == 1:
		year, month, day = config['args'][0].split('/', 2)
		year = int(year)
		month = int(month)
		day = int(day)
		tNow = _MST.localize(datetime(year, month, day))
		tNow = tNow.astimezone(_UTC)
		
	else:
		tNow = datetime.utcnow()
	observer.date = tNow.strftime("%Y/%m/%d %H:%M:%S")
	print "Current time is %s UTC" % tNow.strftime("%Y/%m/%d %H:%M:%S")
	print "CUrrent LST at %s is %s" % (lwa1.name, observer.sidereal_time())
	
	# Load in the sources and compute
	srcs = [ephem.Sun(), ephem.Jupiter()]
	for line in _srcs:
		srcs.append( ephem.readdb(line) )
	for i in xrange(len(srcs)):
		srcs[i].compute(observer)
	
	if not config['positionMode']:
		#
		# Standard prediction output
		#

		# Header
		print ""
		print "%-10s  %-23s  %-23s  %-23s  %-7s" % ("Source", "Next Rise", "Next Transit", "Next Set", "Up Now?")
		print "="*(10+2+23+2+23+2+23+2+7)
	
		# List
		for src in srcs:		
			isUp = False
			if src.alt > 0:
				isUp = True
		
			try:
				nR = str(observer.next_rising(src, tNow.strftime("%Y/%m/%d %H:%M:%S")))
				nR = _UTC.localize( datetime.strptime(nR, "%Y/%m/%d %H:%M:%S") )
				nR = nR.astimezone(_MST)
			except ephem.AlwaysUpError:
				nR = None
		
			nT = str(observer.next_transit(src, start=tNow.strftime("%Y/%m/%d %H:%M:%S")))
			nT = _UTC.localize( datetime.strptime(nT, "%Y/%m/%d %H:%M:%S") )
			nT = nT.astimezone(_MST)
		
			try:
				nS = str(observer.next_setting(src, tNow.strftime("%Y/%m/%d %H:%M:%S")))
				nS = _UTC.localize( datetime.strptime(nS, "%Y/%m/%d %H:%M:%S") )
				nS = nS.astimezone(_MST)
			except ephem.AlwaysUpError:
				nS = None
		
		
			try:
				print "%-10s  %-23s  %-23s  %-23s  %-7s" % (src.name, nR.strftime("%Y/%m/%d %H:%M:%S %Z"), nT.strftime("%Y/%m/%d %H:%M:%S %Z"), nS.strftime("%Y/%m/%d %H:%M:%S %Z"), "*" if isUp else "")
			except AttributeError:
				print "%-10s  %-23s  %-23s  %-23s  %-7s" % (src.name, '---', nT.strftime("%Y/%m/%d %H:%M:%S %Z"), '---', "*" if isUp else "")
	else:
		#
		# Position mode
		#

		# Header
		print ""
		print "%-10s  %-9s  %-9s  %-7s" % ("Source", "  Azimuth", "Elevation", "Rising?")
		print "="*(10+2+9+2+9+2+7)
	
		# List
		for src in srcs:		
			if src.alt <= 0:
				continue

			isRising = False
			if src.az < math.pi:
				isRising = True

			print "%-10s  %9.2f  %9.2f  %7s" % (src.name, src.az*180/math.pi, src.alt*180/math.pi, "Yes" if isRising else "")

if __name__ == "__main__":
	main(sys.argv[1:])
	
