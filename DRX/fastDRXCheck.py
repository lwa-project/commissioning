#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Check the time times in a DRX file for flow in a more intellegent fashion.
This script also allows DRX files to be split at the boundaries of time
flow problems.

$Rev$
$LastChangedBy$
$LastChangedDate$
"""

import os
import sys
import math
import time
import getopt

from lsl.reader import drx, errors


def usage(exitCode=None):
    print """fastDRXCheck.py - Read in a DRX file and identify byte ranges that are
contiguous

Usage: fastDRXCheck.py [OPTIONS] file

Options:
-h, --help             Display this help information
-v, --verbose          Be verbose (default = no)
-l, --loose            Do not require exact time flow, good for LWA-SV files
-m, --min-frames       Minimum number of frames to consider (default = 4096)
-s, --split            Split out sections that are valid (default = no)
-k, --keep             When splitting, only work the N largest sections
                       (default = all)
-b, --brief            Use "brief" filenames (default = long filenames that 
                       contain the section bytes range)
"""
    
    if exitCode is not None:
        sys.exit(exitCode)
    else:
        return True


def parseOptions(args):
    config = {}
    # Command line flags - default values
    config['verbose'] = False
    config['strict'] = True
    config['minFrames'] = 4096
    config['split'] = False
    config['keep'] = -1
    config['brief'] = False
    config['args'] = []
    
    # Read in and process the command line flags
    try:
        opts, args = getopt.getopt(args, "hvlm:sk:b", ["help", "verbose", "loose", "min-frames=", "split", "keep=", "brief"])
    except getopt.GetoptError, err:
        # Print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        usage(exitCode=2)
        
    # Work through opts
    for opt, value in opts:
        if opt in ('-h', '--help'):
            usage(exitCode=0)
        elif opt in ('-v', '--verbose'):
            config['verbose'] = True
        elif opt in ('-l', '--loose'):
            config['strict'] = False
        elif opt in ('-m', '--min-frames'):
            config['minFrames'] = int(value, 10)
        elif opt in ('-s', '--split'):
            config['split'] = True
        elif opt in ('-k', '--keep'):
            config['keep'] = int(value, 10)
        elif opt in ('-b', '--brief'):
            config['brief'] = True
        else:
            assert False
            
    # Add in arguments
    config['args'] = args
    
    # Return configuration
    return config


def identify_section(fh, start=0, stop=-1, strict=True, min_frames=4096, verbose=True):
    if stop <= start:
        stop = os.path.getsize(fh.name)
    fh.seek(start)
    
    # Report
    if verbose:
        print "Working on %i to %i of '%s'..." % (start, stop, os.path.basename(fh.name))
        
    # Make sure this is enough to work with
    if stop-start < drx.FrameSize*min_frames:
        if verbose:
            print "  too small for analysis, skipping"
        return None
        
    # Align on the start of a Mark5C packet...
    while True:
        try:
            junkFrame = drx.readFrame(fh)
            try:
                # ... that has a valid decimation
                srate = junkFrame.getSampleRate()
                break
            except ZeroDivisionError:
                pass
        except errors.syncError:
            fh.seek(-drx.FrameSize+1, 1)
    fh.seek(-drx.FrameSize, 1)
    # ... and save that location
    frame_begin = junkFrame.data.timeTag
    file_begin = fh.tell()
    if verbose:
        print "  start @ %i with %i" % (file_begin, frame_begin)
        
    # Find the last valid Mark5C packet...
    fh.seek(stop-drx.FrameSize)
    while True:
        try:
            junkFrame = drx.readFrame(fh)
            try:
                # ... that has a valid decimation
                srate = junkFrame.getSampleRate()
                break
            except ZeroDivisionError:
                pass
        except errors.syncError:
            fh.seek(-drx.FrameSize-1, 1)
    fh.seek(-drx.FrameSize, 1)
    # ... and save that location
    frame_end = junkFrame.data.timeTag
    file_end = fh.tell() + drx.FrameSize
    if verbose:
        print "  stop  @ %i with %i" % (file_end, frame_end)
        
    # Get how much the timetags should change and other basic information
    fh.seek(file_begin)
    ids = []
    for i in xrange(24*8):
        junkFrame = drx.readFrame(fh)
        b,t,p = junkFrame.parseID()
        id = (t,p)
        if id not in ids:
            ids.append(id)
    ttStep = 4096*junkFrame.header.decimation
    if verbose:
        print "  %i frames with a timetag step of %i" % (len(ids), ttStep)
        
    # Difference
    nBytes = file_end - file_begin
    nFrames = nBytes / drx.FrameSize
    ttDiffFound = frame_end - frame_begin
    ttDiffExpected = nFrames / len(ids) * ttStep
    if verbose:
        print "  -> found timetag difference of    %i" % ttDiffFound
        print "  -> expected timetag difference is %i" % ttDiffExpected
        
    # Decide what to do
    if abs(ttDiffFound - ttDiffExpected) > ttStep*(1-strict):
        if verbose:
            print "  ====> mis-match, subsampling"
        file_middle = file_begin + (nFrames / 2) * drx.FrameSize
        parts0 = identify_section(fh, file_begin, file_middle, strict=strict, min_frames=min_frames, verbose=verbose)
        parts1 = identify_section(fh, file_middle, file_end, strict=strict, min_frames=min_frames, verbose=verbose)
        
    else:
        if verbose:
            print "  ====> good, done"
        parts0 = [[file_begin, file_end],]
        parts1 = None
        
    # Sort and merge
    partsList = []
    for parts in (parts0, parts1):
        if parts is None:
            continue
        for part in parts:
            partsList.append( part )
    partsList.sort()
        
    # Merge
    parts = []
    if len(partsList) > 0:
        parts.append( partsList[0] )
        for part in partsList[1:]:
            if part[0] == parts[-1][1]:
                parts[-1][1] = part[1]
            else:
                parts.append(part)
                
    return parts


def fine_tune_boundary_start(fh, start, max_frames=4, verbose=True):
    fh.seek(start)
    
    # Report
    if verbose:
        print "Tuning boundary at %i of '%s'..." % (start, os.path.basename(fh.name))
        
    # Align on the start of a Mark5C packet...
    while True:
        try:
            junkFrame = drx.readFrame(fh)
            try:
                # ... that has a valid decimation
                srate = junkFrame.getSampleRate()
                break
            except ZeroDivisionError:
                pass
        except errors.syncError:
            fh.seek(-drx.FrameSize+1, 1)
    fh.seek(-drx.FrameSize, 1)
    # ... and save that location
    frame_begin = junkFrame.data.timeTag
    file_begin = fh.tell()
    if verbose:
        print "  start @ %i with %i" % (file_begin, frame_begin)
        
    # Get how much the timetags should change and other basic information
    fh.seek(file_begin)
    ids = []
    for i in xrange(24*8):
        junkFrame = drx.readFrame(fh)
        b,t,p = junkFrame.parseID()
        id = (t,p)
        if id not in ids:
            ids.append(id)
    ttStep = 4096*junkFrame.header.decimation
    if verbose:
        print "  %i frames with a timetag step of %i" % (len(ids), ttStep)
        
    # Load in the times to figure out what to do
    fh.seek(file_begin)
    timetags = []
    for i in xrange(max_frames):
        junkFrame = drx.readFrame(fh)
        timetags.append( junkFrame.data.timeTag )
    skips = [timetags[i]-timetags[i-1] for i in xrange(1, max_frames)]
    try:
        offset = min([skips.index(0), skips.index(ttStep)])
    except ValueError:
        offset = skips.index(0)
    if verbose:
        print "  -> shifting boundary by %i frame(s)" % offset
    start += drx.FrameSize*offset
    return start


def getBestSize(value, powerOfTwo=True):
    scale = 1024. if powerOfTwo else 1000.
    if value >= 0.96*scale**4:
        value = value / scale**4
        unit = 'T'
    elif value >= 0.96*scale**3:
        value = value / scale**3
        unit = 'G'
    elif value >= 0.96*scale**2:
        value = value / scale**2
        unit = 'M'
    elif value >= 0.96*scale**1:
        value = value / scale**1
        unit = 'k'
    else:
        unit = ''
    return value,unit


def main(args):
    # Parse the command line
    config = parseOptions(args)
    filename = config['args'][0]
    
    # Determine how to print out what we find
    scale = math.log10(os.path.getsize(filename))
    scale = int(math.ceil(scale))
    fmt = '  %%%ii, %%%ii -> %%7.3f %%sB or %%7.3f %%sframes' % (scale, scale)
    
    # Figure out the parts
    fh = open(filename, 'rb')
    parts = identify_section(fh, strict=config['strict'], min_frames=config['minFrames'], verbose=config['verbose'])
    if parts is None:
        print "No valid byte ranges found, exiting"
        sys.exit(1)
        
    # Fine tune the boundaries
    for p,part in enumerate(parts):
        start, stop = part
        start = fine_tune_boundary_start(fh, start, verbose=config['verbose'])
        parts[p][0] = start
    
    # Report
    print "Valid Byte Ranges:"
    valid = 0
    rank = []
    for p,part in enumerate(parts):
        start, stop = part
        size = stop - start
        frames = size / drx.FrameSize
        valid += size
        rank.append( [size, start, stop, p] )
        
        s,su = getBestSize(size, powerOfTwo=True)
        f,fu = getBestSize(frames, powerOfTwo=False)
        print fmt % (start, stop, s, su, f, fu)
    print "-> %.1f%% contiguous in %i frame blocks" % (100.0*valid/os.path.getsize(filename), config['minFrames'])
    
    if config['split']:
        print " "
        
        rank.sort(reverse=True)
        if config['keep'] >= 1:
            rank = rank[:config['keep']]
            
        if config['brief']:
            scale = math.log10(len(parts))
            scale = int(math.ceil(scale))
            fmt = '{0:s}-S{3:0%id}' % scale
        else:
            fmt = '{0:s}-{1:0%id}-{2:00%id}' % (scale, scale)
            
        for i,(size,start,stop,section) in enumerate(rank):
            outname = fmt.format(os.path.basename(filename), start, stop, section)
            print "Working on section #%i..." % (i+1)
            print "  Filename: %s" % outname
            
            #size = 4096**2 * 10
            #stop = start + size
            
            t0 = time.time()
            fh.seek(start)
            oh = open(outname, 'wb')
            nBytesRead = size
            for sl in [drx.FrameSize*2**i for i in range(16)[::-1]]:
                while nBytesRead >= sl:
                    temp = fh.read(sl)
                    oh.write(temp)
                    nBytesRead -= sl
            oh.close()
            t1 = time.time()
            s,su = getBestSize(os.path.getsize(outname), powerOfTwo=True)
            print "  Copied %.3f %sB in %.3f s (%.3f MB/s)" % (s, su, t1-t0, os.path.getsize(outname)/1024.0**2/(t1-t0))


if __name__ == "__main__":
    main(sys.argv[1:])
    
