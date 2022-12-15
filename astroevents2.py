#!/usr/bin/env python3

"""
New take on astroevents using PyEphem for calculations.  It can also take a
date in the form of YYYY/MM/DD from the command line to use a as base for 
its calculations.
"""

# Python2 compatibility
from __future__ import print_function, division
try:
    range = xrange
except NameError:
    pass
    
import os
import sys
import pytz
import math
import ephem
import argparse
from datetime import datetime

from lsl.common import stations
from lsl.misc import parser as aph


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
        "CasA,f|J,23:23:27.94,+58:48:42.4,1", ]


def main(args):
    # Set the station
    if args.lwasv:
        station = stations.lwasv
    elif args.ovrolwa:
        station = stations.lwa1
        station.name = 'OVRO-LWA'
        station.lat, station.lon, station.elev = ('37.23977727', '-118.2816667', 1182.89)
    else:
        station = stations.lwa1
    observer = station.get_observer()
    print("Current site is %s at lat %s, lon %s" % (station.name, observer.lat, observer.long))
    
    # Set the current time so we can find the "next" transit.  Go ahead
    # and report the time and the current LST (for reference)
    if args.date is not None and args.time is not None:
        year, month, day = args.date.split('/', 2)
        year = int(year)
        month = int(month)
        day = int(day)
        
        hour, minute, second = args.time.split(':', 2)
        hour = int(hour)
        minute = int(minute)
        second = int(float(second))
        
        tNow = _MST.localize(datetime(year, month, day, hour, minute, second))
        tNow = tNow.astimezone(_UTC)
        
    elif args.date is not None:
        year, month, day = args.date.split('/', 2)
        year = int(year)
        month = int(month)
        day = int(day)
        tNow = _MST.localize(datetime(year, month, day))
        tNow = tNow.astimezone(_UTC)
        
    else:
        tNow = _UTC.localize(datetime.utcnow())
    observer.date = tNow.strftime("%Y/%m/%d %H:%M:%S")
    print("Current time is %s" % tNow.astimezone(_MST).strftime("%Y/%m/%d %H:%M:%S %Z"))
    print("                %s" % tNow.astimezone(_UTC).strftime("%Y/%m/%d %H:%M:%S %Z"))
    print("Current LST at %s is %s" % (station.name, observer.sidereal_time()))
    
    # Load in the sources and compute
    srcs = [ephem.Sun(), ephem.Jupiter()]
    for line in _srcs:
        srcs.append( ephem.readdb(line) )
    for i in range(len(srcs)):
        srcs[i].compute(observer)
        
    if not args.position_mode:
        #
        # Standard prediction output
        #
        
        # Header
        print("")
        print("%-10s  %-23s  %-23s  %-23s  %-7s" % ("Source", "Next Rise", "Next Transit", "Next Set", "Up Now?"))
        print("="*(10+2+23+2+23+2+23+2+7))
        
        # List
        for src in srcs:		
            isUp = False
            if src.alt > 0:
                isUp = True
                
            try:
                nR = str(observer.next_rising(src, tNow.strftime("%Y/%m/%d %H:%M:%S")))
                nR = _UTC.localize( datetime.strptime(nR, "%Y/%m/%d %H:%M:%S") )
                if not args.utc:
                    nR = nR.astimezone(_MST)
            except ephem.AlwaysUpError:
                nR = None
                
            nT = str(observer.next_transit(src, start=tNow.strftime("%Y/%m/%d %H:%M:%S")))
            nT = _UTC.localize( datetime.strptime(nT, "%Y/%m/%d %H:%M:%S") )
            if not args.utc:
                nT = nT.astimezone(_MST)
                
            try:
                nS = str(observer.next_setting(src, tNow.strftime("%Y/%m/%d %H:%M:%S")))
                nS = _UTC.localize( datetime.strptime(nS, "%Y/%m/%d %H:%M:%S") )
                if not args.utc:
                    nS = nS.astimezone(_MST)
            except ephem.AlwaysUpError:
                nS = None
                
                
            try:
                print("%-10s  %-23s  %-23s  %-23s  %-7s" % (src.name, nR.strftime("%Y/%m/%d %H:%M:%S %Z"), nT.strftime("%Y/%m/%d %H:%M:%S %Z"), nS.strftime("%Y/%m/%d %H:%M:%S %Z"), "*" if isUp else ""))
            except AttributeError:
                print("%-10s  %-23s  %-23s  %-23s  %-7s" % (src.name, '---', nT.strftime("%Y/%m/%d %H:%M:%S %Z"), '---', "*" if isUp else ""))
    else:
        #
        # Position mode
        #
        
        # Header
        print("")
        print("%-10s  %-9s  %-9s  %-7s" % ("Source", "  Azimuth", "Elevation", "Rising?"))
        print("="*(10+2+9+2+9+2+7))
        
        # List
        for src in srcs:		
            if src.alt <= 0:
                continue
                
            isRising = False
            if src.az < math.pi:
                isRising = True
                
            print("%-10s  %9.2f  %9.2f  %7s" % (src.name, src.az*180/math.pi, src.alt*180/math.pi, "Yes" if isRising else ""))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='new take on the astroevents.py script included in LSL that can make calculations for a different day',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
    parser.add_argument('date', type=aph.date, nargs='?',
                        help='MST/MDT date for calculation in YYYY/MM/DD')
    parser.add_argument('time', type=aph.time, nargs='?',
                        help='MST/MDT time for calculation in HH:MM:SS')
    parser.add_argument('-u', '--utc', action='store_true',
                        help='display rise, transit, and set times in UTC instead of MST/MDT')
    parser.add_argument('-p', '--position-mode', action='store_true',
                        help='display the azimuth and elevation of sources above the horizon')
    sgroup = parser.add_mutually_exclusive_group(required=False)
    sgroup.add_argument('-s', '--lwasv', action='store_true',
                        help='compute for LWA-SV instead of LWA1')
    sgroup.add_argument('-o', '--ovrolwa', action='store_true',
                        help='compute for OVRO-LWA instead of LWA1')
    args = parser.parse_args()
    main(args)
    
