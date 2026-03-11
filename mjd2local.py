#!/usr/bin/env python3

"""
For a given MJD value or list of MJD values, return the range of local times 
associated with that MJD.
"""

import sys
import argparse
from datetime import datetime, timezone
try:
    import zoneinfo
except ImportError:
    from backports import zoneinfo
    
from lsl.common.mcs import mjdmpm_to_datetime
from lsl.misc import parser as aph

_MST = zoneinfo.ZoneInfo('America/Denver')


def main(args):
    otz = _MST
    if args.utc:
        otz = timezone.utc
        
    if not args.pairs:
        for arg in args.mjd:
            mjd1 = int(arg)
            mjd2 = float(mjd1) + 0.99999

            d1 = mjdmpm_to_datetime(mjd1, 0, tz=timezone.utc)
            d1  = d1.astimezone(otz)

            d2 = mjdmpm_to_datetime(mjd2, 0, tz=timezone.utc)
            d2  = d2.astimezone(otz)
            
            tzname = d1.strftime('%Z')
            
            print("MJD: %i" % mjd1)
            print("%s: %s to %s" % (tzname, d1.strftime("%B %d, %Y at %H:%M:%S %Z"), d2.strftime("%B %d, %Y at %H:%M:%S %Z")))
    else:
        for arg in zip(args.mjd[0::2], args.mjd[1::2]):
            mjd, mpm = [int(i) for i in arg]
            d = mjdmpm_to_datetime(mjd, mpm, tz=timezone.utc)
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
    
