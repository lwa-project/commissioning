#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
splitObservations.py - Read in a DRX/HDF5 watefall file and split out various 
observations.  The observations can be split by:
* Target Source
* Observation ID number
or the script can be used to just list the observations within an HDF5 file.

Usage:
./splitObservations.py [OPTIONS] file
"""

import os
import sys
import h5py
import numpy
import argparse
from datetime import datetime


def main(args):
    h = h5py.File(args.filename, 'r')
    
    if args.list:
        for obsName in h.keys():
            obs = h.get(obsName, None)
            
            # Load in the information we need to calculate the pseudo-spectral kurtosis
            target = obs.attrs['TargetName']
            mode = obs.attrs['TrackingMode']
            tInt = obs.attrs['tInt']
            LFFT = obs.attrs['LFFT']
            srate = obs.attrs['sample_rate']
            
            print "Observation #%i" % int(obsName.replace('Observation', ''))
            print "  Target: %s" % target
            print "  Mode: %s" % mode
            print "  Sample Rate: %.1f Hz" % srate
            print "  LFFT: %i" % LFFT
            print "  tInt: %.3f s" % tInt
            
    else:
        if args.source:
            # Create a list of sources and which observation ID they correspond to
            sources = {}
            for obsName in h.keys():
                obsID = int(obsName.replace('Observation', ''))
                obs = h.get(obsName, None)
                
                try:
                    sources[obs.attrs['TargetName']].append(obsID)
                except KeyError:
                    sources[obs.attrs['TargetName']] = [obsID,]
                    
            # Loop over those sources and create a new HDF5 file for each
            for source,obsIDs in sources.iteritems():
                outname = os.path.split(args.filename)[1]
                outname = os.path.splitext(outname)[0]
                outname = "%s-%s.hdf5" % (outname, source.replace(' ', '').replace('/','').replace('&','and'))
                
                if os.path.exists(outname):
                    if not args.force:
                        yn = raw_input("WARNING: '%s' exists, overwrite? [Y/n] " % outname)
                    else:
                        yn = 'y'
                        
                    if yn not in ('n', 'N'):
                        os.unlink(outname)
                    else:
                        print "WARNING: output file '%s' already exists, skipping" % outname
                        continue
                        
                hOut = h5py.File(outname, mode='a')
                for name in h.attrs.keys():
                    hOut.attrs[name] = h.attrs[name]
                for i,obsID in enumerate(obsIDs):
                    h.copy("Observation%i" % obsID, hOut, name='Observation%i' % (i+1))
                hOut.attrs['FileCreation'] = datetime.utcnow().strftime("UTC %Y/%m/%d %H:%M:%S")
                hOut.close()
                
        else:
            # Loop over all of the observations and create a new HDF5 file for each one
            for obsName in h.keys():
                obsID = int(obsName.replace('Observation', ''))
                
                outname = os.path.split(args.filename)[1]
                outname = os.path.splitext(outname)[0]
                outname = "%s-%i.hdf5" % (outname, obsID)
                
                if os.path.exists(outname):
                    if not args.force:
                        yn = raw_input("WARNING: '%s' exists, overwrite? [Y/n] " % outname)
                    else:
                        yn = 'y'
                        
                    if yn not in ('n', 'N'):
                        os.unlink(outname)
                    else:
                        print "WARNING: output file '%s' already exists, skipping" % outname
                        continue
                        
                hOut = h5py.File(outname, 'a')
                for name in h.attrs.keys():
                    hOut.attrs[name] = h.attrs[name]
                h.copy(obsName, hOut, name='Observation1')
                hOut.attrs['FileCreation'] = datetime.utcnow().strftime("UTC %Y/%m/%d %H:%M:%S")
                hOut.close()
                
    h.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='read in DRX/HDF5 waterfall file split out various observations', 
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
    parser.add_argument('filename', type=str, 
                        help='filename to decimate')
    parser.add_argument('-l', '--list', action='store_true', 
                        help='list source names')
    parser.add_argument('-s', '--source', action='store_true', 
                        help='split by source name instead of observation ID')
    parser.add_argument('-f', '--force', action='store_true', 
                        help='force overwritting of existing HDF5 files')
    args = parser.parse_args()
    main(args)
    
