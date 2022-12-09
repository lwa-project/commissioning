#!/usr/bin/env python3

# Python3 compatiability
from __future__ import print_function, division
import sys
if sys.version_info > (3,):
    xrange = range
    
import os
import sys
import pytz
import math
import ephem
import getopt
from datetime import datetime, timedelta

from lsl.common.stations import lwa1
from lsl.common.mcs import datetime_to_mjdmpm
from lsl.reader.tbn import FRAME_SIZE as tbnFRAME_SIZE
from lsl.reader.tbn import FILTER_CODES as tbn_filters


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
         "J2145-0750,f|J,21:45:50.47,-07:50:18.3,1",
         "J2339-0533,f|J,23:39:38.75,-05:33:05.3,1"]


def usage(exitCode=None):
    print("""tbnCalibRecordCommand.py - Get a DR REC command to record 100 kS/s 
TBN (filter #7) for two hours centered on the transit of Cyg A.

Note:  Input dates are in UTC.

Usage: tbnCalibRecordCommand.py [OPTIONS] [YYYY/MM/DD]

Options:
-h, --help             Display this help information
-s, --source           Source to center on (default = CygA)
-d, --duration         Duration of recording in seconds 
                    (default = 7200)
-f, --filter           TBN filter (default = 7)
""")
    
    if exitCode is not None:
        sys.exit(exitCode)
    else:
        return True


def parseOptions(args):
    config = {}
    config['source'] = 'CygA'
    config['duration'] = 7200
    config['filter'] = 7
    
    # Read in and process the command line flags
    try:
        opts, args = getopt.getopt(args, "hs:d:f:", ["help", "source=", "duration=", "filter="])
    except getopt.GetoptError as err:
        # Print help information and exit:
        print(str(err)) # will print something like "option -a not recognized"
        usage(exitCode=2)
        
    # Work through opts
    for opt, value in opts:
        if opt in ('-h', '--help'):
            usage(exitCode=0)
        elif opt in ('-s', '--source'):
            config['source'] = value
        elif opt in ('-d', '--duration'):
            config['duration'] = float(value)
        elif opt in ('-f', '--filter'):
            config['filter'] = int(value)
        else:
            assert False
            
    # Validate the filter
    if config['filter'] < 1 or config['filter'] > 7:
        raise ValueError("Invalid TBN filter code: %i" % config['filter'])
        
    # Add in arguments
    config['args'] = args
    
    # Return configuration
    return config


def main(args):
    # Parse the command line
    config = parseOptions(args)
    
    # Get LWA-1
    observer = lwa1.get_observer()
    print("Current site is %s at lat %s, lon %s" % (lwa1.name, observer.lat, observer.long))
    
    # Set the current time so we can find the "next" transit.  Go ahead
    # and report the time and the current LST (for reference)
    if len(config['args']) == 1:
        config['args'][0] = config['args'][0].replace('-', '/')
        year, month, day = config['args'][0].split('/', 2)
        year = int(year)
        month = int(month)
        day = int(day)
        tNow = _UTC.localize(datetime(year, month, day))
        
    else:
        tNow = _UTC.localize(datetime.utcnow())
    observer.date = tNow.strftime("%Y/%m/%d %H:%M:%S")
    print("Current time is %s" % tNow.astimezone(_UTC).strftime("%Y/%m/%d %H:%M:%S %Z"))
    print("Current LST at %s is %s" % (lwa1.name, observer.sidereal_time()))
    
    # Load in the sources and compute
    srcs = [ephem.Sun(), ephem.Jupiter()]
    for line in _srcs:
        srcs.append( ephem.readdb(line) )
    for i in xrange(len(srcs)):
        srcs[i].compute(observer)
        
    #
    # Standard prediction output
    #
    
    # Header
    print("")
    print("%-10s  %-23s" % ("Source", "Next Transit", ))
    print("="*(10+2+23))
    
    # List
    found = False
    for src in srcs:
        if src.name.lower() == config['source'].lower():
            found = True
            
            nT = str(observer.next_transit(src, start=tNow.strftime("%Y/%m/%d %H:%M:%S")))
            nT = _UTC.localize( datetime.strptime(nT, "%Y/%m/%d %H:%M:%S") )
            
            print("%-10s %-23s" % (src.name, nT.strftime("%Y/%m/%d %H:%M:%S %Z")))
            print("%-10s %-23s" % ("", nT.astimezone(_MST).strftime("%Y/%m/%d %H:%M:%S %Z")))
            print(" ")
            break
            
    if found:
        startRec = nT - timedelta(seconds=int(round(config['duration']/2.0)))
        mjd, mpm = datetime_to_mjdmpm(startRec)
        dur = int(round(config['duration']*1000))
        cmd = 'DR5 REC "%i %i %i TBN_FILT_%i"' % (mjd, mpm, dur, config['filter'])
        
        antpols = 520
        sample_rate = tbn_filters[config['filter']]
        dataRate = 1.0*sample_rate/512*tbnFRAME_SIZE*antpols
        
        print("Recording:")
        print(" Start: %s" % startRec.strftime("%Y/%m/%d %H:%M:%S %Z"))
        print("        %s" % startRec.astimezone(_MST).strftime("%Y/%m/%d %H:%M:%S %Z"))
        print(" Duration: %.3f s" % (dur/1000.0,))
        print(" TBN Filter Code: %i" % config['filter'])
        print(" Data rate: %.2f MB/s" % (dataRate/1024**2,))
        print(" Data volume: %.2f GB" % (dataRate*dur/1000.0/1024**3,))
        print(" ")
        
        print("Data Recorder Command:")
        print(" %s" % cmd)
    else:
        raise RuntimeError("Unknown source '%s'" % config['source'])


if __name__ == "__main__":
    main(sys.argv[1:])
    