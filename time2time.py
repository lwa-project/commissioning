#!/usr/bin/env python3

"""
Convert a local date/time string in the format of "YYYY/MM/DD HH:MM:SS[.SSS]" into 
MJD and MPM UTC values.  If no date/time string is supplied, the current local 
date/time is used.
"""

# Python3 compatiability
from __future__ import print_function, division
import sys
if sys.version_info > (3,):
    xrange = range
    
import sys
import math
import pytz
import argparse
import datetime

from lsl.common.stations import lwa1
from lsl.common.mcs import datetime_to_mjdmpm
from lsl.astro import date as astroDate, get_date as astroGetDate
from lsl.misc import parser as aph


def _getEquinoxEquation(jd):
    """
    Compute the equation of the equinoxes (nutation in right ascension) in 
    hours for the specified Julian Date.
    
    From:
    http://aa.usno.navy.mil/faq/docs/GAST.php
    """
    
    # Get the number of days since January 1, 2000 @ 12:00 UT
    D = jd - 2451545.0

    # Compute the obliquity
    epsilon = 23.4393 - 0.0000004*D
    
    # Compute the mean longitude of the Sun
    L = 280.47 + 0.98565*D
    
    # Compute the longitude of the Moon's ascending node
    Omega = 125.04 - 0.052954*D
    
    # The nutation in the longitude (hours)
    deltaPsi = -0.000319*math.sin(Omega*math.pi/180) - 0.000024*math.sin(2*L*math.pi/180)

    # Return the equation of the equinoxes
    return deltaPsi * math.cos(epsilon*math.pi/180.0)


def main(args):
    obs = lwa1.get_observer()
    MST = pytz.timezone('US/Mountain')
    UTC = pytz.utc
    
    if args.date is None or args.time is None:
        dt = datetime.datetime.utcnow()
        dt = UTC.localize(dt)
    else:
        year, month, day = args.date.split('/', 2)
        hour, minute, second = args.time.split(':', 2)
        iSeconds = int(float(second))
        mSeconds = int(round((float(second) - iSeconds)*1000000))
        
        if args.utc:
            # UTC
            dt = UTC.localize(datetime.datetime(int(year), int(month), int(day), int(hour), int(minute), iSeconds, mSeconds))
            
        elif args.sidereal:
            # LST
            dt = astroDate(int(year), int(month), int(day), 0, 0, 0)
            jd = dt.to_jd()
            
            ## Get the LST in hours
            LST = int(hour) + int(minute)/60.0 + (iSeconds + mSeconds/1e6)/3600.0
            
            ## Get the Greenwich apparent ST for LST using the longitude of 
            ## the site.  The site longitude is stored as radians, so convert
            ## to hours first.
            GAST = LST - obs.long*12/math.pi
            
            ## Get the Greenwich mean ST by removing the equation of the 
            ## equinoxes (or some approximation thereof)
            GMST = GAST - _getEquinoxEquation(jd)
            
            ## Get the value of D0, days since January 1, 2000 @ 12:00 UT, 
            ## and T, the number of centuries since the year 2000.  The value
            ## of T isn't terribly important but it is nice to include
            D0 = jd - 2451545.0
            T = D0 / 36525.0
            
            ## Solve for the UT hour for this LST and map onto 0 -> 24 hours
            ## From: http://aa.usno.navy.mil/faq/docs/GAST.php
            H  = GMST - 6.697374558 - 0.06570982441908*D0 - 0.000026*T**2
            H /= 1.002737909350795
            while H < 0:
                H += 24/1.002737909350795
            while H > 24:
                H -= 24/1.002737909350795
                
            ## Get the full Julian Day that this corresponds to
            jd += H/24.0
            
            ## Convert the JD back to a time and extract the relevant 
            ## quantities needed to build a datetime instance
            dt = astroGetDate(jd)
            year = dt.years
            month = dt.months
            day = dt.days
            hour = dt.hours
            minute = dt.minutes
            second = int(dt.seconds)
            microsecond = int((dt.seconds - second)*1e6)
            ## Trim the microsecond down to the millisecond level
            microsecond = int(int(microsecond/1000.0)*1000)
            
            ## Localize as the appropriate time zone
            dt = UTC.localize(datetime.datetime(year, month, day, hour, minute, second, microsecond))
            
        else:
            # Mountain time
            dt = MST.localize(datetime.datetime(int(year), int(month), int(day), int(hour), int(minute), iSeconds, mSeconds))
            
        dt = dt.astimezone(UTC)
        
    obs.date = dt.astimezone(UTC).strftime("%Y/%m/%d %H:%M:%S.%f")
    mjd, mpm = datetime_to_mjdmpm(dt)
    
    print("Localtime: %s" % dt.astimezone(MST).strftime("%B %d, %Y at %H:%M:%S %Z"))
    print("UTC: %s" % dt.astimezone(UTC).strftime("%B %d, %Y at %H:%M:%S %Z"))
    print("LST: %s" % obs.sidereal_time())
    print("MJD: %i" % mjd)
    print("MPM: %i" % mpm)
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='convert a local date/time string into MJD and MPM UTC values',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
    parser.add_argument('date', type=aph.date, nargs='?',
                        help='local date in YYYY/MM/DD')
    parser.add_argument('time', type=aph.time, nargs='?',
                        help='local time in HH:MM:SS[.SSS]')
    group = parser.add_mutually_exclusive_group(required=False)
    group.add_argument('-s', '--sidereal', action='store_true',
                       help='iput time is in LST, not local')
    group.add_argument('-u', '--utc', action='store_true',
                       help='input time is in UTC, not local')
    args = parser.parse_args()
    main(args)
    