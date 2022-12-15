#!/usr/bin/env python3

"""
For a given MJD value or list of MJD values, return the range of local times 
associated with that MJD.
"""

# Python2 compatiability
from __future__ import print_function, division

import sys
import pytz
import argparse
from datetime import datetime

from lsl.common.mcs import mjdmpm_to_datetime
from lsl.misc import parser as aph

MST = pytz.timezone('US/Mountain')
UTC = pytz.utc


def main(args):
    otz = MST
    if args.utc:
        otz = UTC
        
    if not args.pairs:
        for arg in args.mjd:
            mjd1 = int(arg)
            mjd2 = float(mjd1) + 0.99999

            d1 = mjdmpm_to_datetime(mjd1, 0)
            d1 = UTC.localize(d1)
            d1  = d1.astimezone(otz)

            d2 = mjdmpm_to_datetime(mjd2, 0)
            d2 = UTC.localize(d2)
            d2  = d2.astimezone(otz)
            
            tzname = d1.strftime('%Z')
            
            print("MJD: %i" % mjd1)
            print("%s: %s to %s" % (tzname, d1.strftime("%B %d, %Y at %H:%M:%S %Z"), d2.strftime("%B %d, %Y at %H:%M:%S %Z")))
    else:
        for arg in zip(args.mjd[0::2], args.mjd[1::2]):
            mjd, mpm = [int(i) for i in arg]
            d = mjdmpm_to_datetime(mjd, mpm)
            d = UTC.localize(d)
            d = d.astimezone(otz)
            
            tzname = d.strftime('%Z')
            
            print("MJD: %i, MPM: %i" % (mjd, mpm))
            print("%s: %s" % (tzname, d.strftime("%B %d, %Y at %H:%M:%S %Z")))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='for a given MJD value or list of MJD values, return the range of local times associated with that MJD',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
    parser.add_argument('mjd', type=aph.positive_int, nargs='+',
                        help='local date in YYYY/MM/DD')
    parser.add_argument('-u', '--utc', action='store_true',
                       help='report times in UTC rather than local')
    parser.add_argument('-p', '--pairs', action='store_true',
                        help='interpret the input as MJD, MPM pairs')
    args = parser.parse_args()
    main(args)
    
