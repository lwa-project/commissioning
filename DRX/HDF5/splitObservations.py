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

$Rev$
$LastChangedBy$
$LastChangedDate$
"""

import os
import sys
import h5py
import numpy
import getopt
from datetime import datetime


def usage(exitCode=None):
	print """splitObservations.py - Read in DRX/HDF5 waterfall file split out various 
observations.

Usage: splitObservations.py [OPTIONS] file

Options:
-h, --help                Display this help information
-l, --list                List source names
-s, --source              Split by source name instead of observation 
                          ID
-f, --force               Force overwritting of existing HDF5 files
"""

	if exitCode is not None:
		sys.exit(exitCode)
	else:
		return True


def parseOptions(args):
	config = {}
	# Command line flags - default values
	config['list'] = False
	config['source'] = False
	config['force'] = False
	config['args'] = []

	# Read in and process the command line flags
	try:
		opts, args = getopt.getopt(args, "hlsf", ["help", "list", "source", "force"])
	except getopt.GetoptError, err:
		# Print help information and exit:
		print str(err) # will print something like "option -a not recognized"
		usage(exitCode=2)
	
	# Work through opts
	for opt, value in opts:
		if opt in ('-h', '--help'):
			usage(exitCode=0)
		elif opt in ('-l', '--list'):
			config['list'] = True
		elif opt in ('-s', '--source'):
			config['source'] = True
		elif opt in ('-f', '--force'):
			config['force'] = True
		else:
			assert False
	
	# Add in arguments
	config['args'] = args

	# Return configuration
	return config

def main(args):
	config = parseOptions(args)
	
	filename = config['args'][0]
	h = h5py.File(filename, 'r')
	
	if config['list']:
		for obsName in h.keys():
			obs = h.get(obsName, None)
			
			# Load in the information we need to calculate the pseudo-spectral kurtosis
			target = obs.attrs['TargetName']
			mode = obs.attrs['TrackingMode']
			tInt = obs.attrs['tInt']
			LFFT = obs.attrs['LFFT']
			srate = obs.attrs['sampleRate']
			
			print "Observation #%i" % int(obsName.replace('Observation', ''))
			print "  Target: %s" % target
			print "  Mode: %s" % mode
			print "  Sample Rate: %.1f Hz" % srate
			print "  LFFT: %i" % LFFT
			print "  tInt: %.3f s" % tInt
			
	else:
		if config['source']:
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
				outname = os.path.split(filename)[1]
				outname = os.path.splitext(outname)[0]
				outname = "%s-%s.hdf5" % (outname, source.replace(' ', '').replace('/','').replace('&','and'))
				
				if os.path.exists(outname):
					if not config['force']:
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
				
				outname = os.path.split(filename)[1]
				outname = os.path.splitext(outname)[0]
				outname = "%s-%i.hdf5" % (outname, obsID)
				
				if os.path.exists(outname):
					if not config['force']:
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
	main(sys.argv[1:])
	
