#!/usr/bin/env python

"""
Export select stands from a TBW file to HDF5.

Usage:
./tbw2hdf.py <TBW_filename> <stand_ID> [<stand_ID> [...]]
"""

# Python3 compatiability
from __future__ import print_function, division
import sys
if sys.version_info > (3,):
    xrange = range
    
import os
import sys
import h5py
import ephem
import numpy

from lsl.common.stations import lwa1
from lsl.reader import tbw, tbn
from lsl.reader import errors
from lsl.astro import unix_to_utcjd, DJD_OFFSET


def main(args):
    filename = args[0]
    stands = [int(i) for i in args[1:]]

    antennas = lwa1.antennas
    
    fh = open(filename, "rb")
    nFrames = os.path.getsize(filename) // tbw.FRAME_SIZE
    dataBits = tbw.get_data_bits(fh)
    # The number of ant/pols in the file is hard coded because I cannot figure out 
    # a way to get this number in a systematic fashion
    antpols = len(antennas)
    if dataBits == 12:
        nSamples = 400
    else:
        nSamples = 1200

    # Read in the first frame and get the date/time of the first sample 
    # of the frame.  This is needed to get the list of stands.
    junkFrame = tbw.read_frame(fh)
    fh.seek(0)
    beginTime = sum(junkFrame.time, 0.0)
    beginDate = ephem.Date(unix_to_utcjd(sum(junkFrame.time, 0.0)) - DJD_OFFSET)
    
    # Figure out which digitizers to keep
    toKeep = []
    for a in antennas:
        if a.stand.id in stands:
            toKeep.append( a.digitizer )

    # File summary
    print("Filename: %s" % filename)
    print("Date of First Frame: %s" % str(beginDate))
    print("Ant/Pols: %i" % antpols)
    print("Sample Length: %i-bit" % dataBits)
    print("Frames: %i" % nFrames)
    print("===")
    print("Keeping Stands:")
    for a in toKeep:
        print(" Stand #%3i, pol %i (digitizer %3i)" % (antennas[a-1].stand.id, antennas[a-1].pol, antennas[a-1].digitizer))

    # Skip over any non-TBW frames at the beginning of the file
    i = 0
    junkFrame = tbw.read_frame(fh)
    while not junkFrame.header.is_tbw:
        try:
            junkFrame = tbw.read_frame(fh)
        except errors.SyncError:
            fh.seek(0)
            while True:
                try:
                    junkFrame = tbn.read_frame(fh)
                    i += 1
                except errors.SyncError:
                    break
            fh.seek(-2*tbn.FRAME_SIZE, 1)
            junkFrame = tbw.read_frame(fh)
        i += 1
    fh.seek(-tbw.FRAME_SIZE, 1)
    print("Skipped %i non-TBW frames at the beginning of the file" % i)
    
    # Create the HDF5 file
    outname = os.path.splitext(filename)[0]
    outname = "%shdf5" % outname
    
    f = h5py.File(outname, 'w')
    f.attrs['filename'] = filename
    f.attrs['mode'] = 'TBW'
    f.attrs['station'] = 'LWA-1'
    f.attrs['dataBits'] = dataBits
    f.attrs['startTime'] = beginTime
    f.attrs['startTime_units'] = 's'
    f.attrs['startTime_sys'] = 'unix'
    f.attrs['sample_rate'] = 196e6
    f.attrs['sample_rate_units'] = 'Hz'
    
    ## Create the digitizer to dataset lookup table and the 
    standLookup = {}
    standData = []
    i = 0
    for a in toKeep:
        if a % 2 == 0:
            continue
        
        s = antennas[a-1].stand.id
        standLookup[a] = i
        
        temp = f.create_group('Stand%03i' % s)
        ### Combined status code
        temp.attrs['statusCode'] = antennas[a-1].combined_status
        
        ### Antenna number
        temp.attrs['antennaID'] = antennas[a-1].id
        
        ### Cable information
        temp.attrs['cableID'] = antennas[a-1].cable.id
        temp.attrs['cableLength'] = antennas[a-1].cable.length
        temp.attrs['cableLength_units'] = 'm'
        
        ### Stand location information
        temp.attrs['posX'] = antennas[a-1].stand.x
        temp.attrs['posX_units'] = 'm'
        temp.attrs['posY'] = antennas[a-1].stand.y
        temp.attrs['posY_units'] = 'm'
        temp.attrs['posZ'] = antennas[a-1].stand.z
        temp.attrs['posZ_units'] = 'm'
        
        ### Time series data sets
        temp.attrs['axis0'] = 'time'
        xpol = temp.create_dataset('X', (12000000,), 'i2', chunks=True)
        ypol = temp.create_dataset('Y', (12000000,), 'i2', chunks=True)
        
        standData.append( (xpol, ypol, temp) )
        i += 1
    
    # Go!
    while True:
        # Read in the next frame and anticipate any problems that could occur
        try:
            cFrame = tbw.read_frame(fh)
        except errors.EOFError:
            break
        except errors.SyncError:
            print("WARNING: Mark 5C sync error on frame #%i" % (int(fh.tell())/tbw.FRAME_SIZE-1))
            continue
        if not cFrame.header.is_tbw:
            continue
        
        # Get the DP "stand" ID and the digitizer number
        stand = cFrame.header.id
        aStand = 2*(stand-1)
        digitizer = aStand + 1
        
        # If we don't need it, skip it
        if digitizer not in toKeep:
            continue
        
        # Actually load the data.
        ## Frame count
        count = cFrame.header.frame_count - 1
        ## Which data set
        dataset = standLookup[digitizer]
        ## Load
        standData[dataset][0][count*nSamples:(count+1)*nSamples] = cFrame.payload.data[0,:]
        standData[dataset][1][count*nSamples:(count+1)*nSamples] = cFrame.payload.data[1,:]
        
    fh.close()
    f.close()


if __name__ == "__main__":
    main(sys.argv[1:])
    