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
import ephem
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
         "B0950+08,f|J,09:53:09.31,+07:55:35.8,1",
         "B1133+16,f|J,11:36:03.25,+15:51:04.5,1",
         "B1919+21,f|J,19:21:44.80,+21:53:01.8,1",]


def main(args):	
	# Get LWA-1
	observer = lwa1.getObserver()
	print "Current site is %s at lat %s, lon %s" % (lwa1.name, observer.lat, observer.long)
	
	# Set the current time so we can find the "next" transit.  Go ahead
	# and report the time and the current LST (for reference)
	if len(args) == 1:
		year, month, day = args[0].split('/', 2)
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
	
	# Header
	print ""
	print "%-8s  %-23s  %-23s  %-23s  %-7s" % ("Source", "Next Rise", "Next Transit", "Next Set", "Up Now?")
	print "="*(8+2+23+2+23+2+23+2+7)
	
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
			print "%-8s  %-23s  %-23s  %-23s  %-7s" % (src.name, nR.strftime("%Y/%m/%d %H:%M:%S %Z"), nT.strftime("%Y/%m/%d %H:%M:%S %Z"), nS.strftime("%Y/%m/%d %H:%M:%S %Z"), "*" if isUp else "")
		except AttributeError:
			print "%-8s  %-23s  %-23s  %-23s  %-7s" % (src.name, '---', nT.strftime("%Y/%m/%d %H:%M:%S %Z"), '---', "*" if isUp else "")


if __name__ == "__main__":
	main(sys.argv[1:])
	