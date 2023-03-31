#!/usr/bin/env python3

"""
Create an eye diagram for some portion of a TBN or DRX file.
"""

import os
import sys
import numpy
import argparse

from lsl.common.dp import fS
from lsl.reader import errors, tbn, drx
from lsl.misc import parser as aph

from matplotlib import pyplot as plt


def main(args):
    # Filename
    filename = args.filename
    
    # The observation and sinusoid frequencies as well as the plotting time
    cFreq = args.freq
    iFreq = args.input_freq
    tInt = args.time
    
    # Reader
    rdr = tbn
    
    #
    # Step 1:  Read in the file and process it
    #
    
    # File size
    sizeB = os.path.getsize(filename)
    
    # Open the file
    fh = open(filename, 'rb')

    # Get the sampel rate and number of stands for each pol
    sample_rate = rdr.get_sample_rate(fh)
    nFrames = numpy.array(rdr.get_frames_per_obs(fh))
    nCaptures = sizeB // rdr.FRAME_SIZE // nFrames.sum()

    print("Filename:    %s" % filename)
    print("Size:        %.1f MB" % (float(sizeB)/1024/1024))
    print("Captures:    %i" % nCaptures)
    print("Sample Rate: %.2f kHz" % (sample_rate/1000.0))
    print("===")

    frame = rdr.read_frame(fh)
    fh.seek(0)
    
    # Offset in frames for beampols beam/tuning/pol. sets
    offset = int(args.skip * sample_rate / frame.payload.data.size * nFrames.sum())
    offset = int(1.0 * offset / nFrames.sum()) * nFrames.sum()
    args.skip = 1.0 * offset / nFrames.sum() * frame.payload.data.size / sample_rate
    fh.seek(offset*rdr.FRAME_SIZE)
    
    nCaptures = (sizeB - offset*rdr.FRAME_SIZE) / rdr.FRAME_SIZE / nFrames.sum()
    
    # Compute the integration time and the number of frames per stand per 
    # integration
    tInt = int(round( tInt * sample_rate)) / sample_rate
    fInt = int(tInt * sample_rate / frame.payload.data.size)
    
    # More output
    print("Keeping only:", args.keep)
    print("Skipping: %.3f s" % args.skip)
    print("Integration Time: %.3f s" % tInt)
    print("Number of integrations in file: %i" % (nCaptures/fInt))
    print(" ")
    
    # Go...
    dtime = numpy.zeros((nFrames.sum(), fInt*frame.payload.data.size), dtype=numpy.float64)
    data = numpy.zeros((nFrames.sum(), fInt*frame.payload.data.size), dtype=numpy.complex64)
    
    count = {}
    standMapper = []
    masterCount = 0
    tStart = 0
    data *= 0
    for f in range(fInt*nFrames.sum()):
        # Read in the next frame and anticipate any problems that could occur
        try:
            cFrame = rdr.read_frame(fh, verbose=False)
        except errors.EOFError:
            break
        except errors.SyncError:
            #print("WARNING: Mark 5C sync error on frame #%i" % (int(fh.tell())/rdr.FRAME_SIZE-1))
            continue
        
        if f == 0:
            tStart = cFrame.time

        try:
            beam,tune,pol = cFrame.id
            aStand = 4*(beam-1) + 2*(tune-1) + pol
        except:
            stand,pol = cFrame.header.id
            aStand = 2*(stand-1)+pol
        
        if aStand not in standMapper:
            standMapper.append(aStand)
            oStand = 1*aStand
            aStand = standMapper.index(aStand)
            try:
                print("Mapping stand %i, pol. %1i (%2i) to array index %3i" % (stand, pol, oStand, aStand))
            except:
                print("Mapping beam %i, tune. %1i, pol. %1i (%2i) to array index %3i" % (beam, tune, pol, oStand, aStand))
        else:
            aStand = standMapper.index(aStand)
        
        if aStand not in count.keys():
            count[aStand] = 0

        dtime[aStand, count[aStand]*cFrame.payload.data.size:(count[aStand]+1)*cFrame.payload.data.size] = 512*count[aStand]/sample_rate + 1.0 / sample_rate * numpy.arange(0.0, cFrame.payload.data.size, dtype=numpy.float64)
        data[aStand, count[aStand]*cFrame.payload.data.size:(count[aStand]+1)*cFrame.payload.data.size] = cFrame.payload.data
        
        count[aStand] = count[aStand] + 1
        masterCount = masterCount + 1
            
    #
    # Step 2: Create the eye diagram
    #
    period = 1.0 / numpy.abs(iFreq - cFreq)
    dtime = (dtime - dtime.min()) % period / period
    
    endPt = data.shape[1]/8
    print(endPt / sample_rate / period)
    
    fig = plt.figure()
    if args.polar:
        sortedMapper = sorted(standMapper)
        for k, aStand in enumerate(sortedMapper):
            s = standMapper.index(aStand)
            if s+1 not in args.keep:
                continue
            
            ax1 = fig.add_subplot(2, 2, k+1, polar=True)
            
            ax1.plot(dtime[s,:]*2*numpy.pi, data[s,:].real - data[s,:].real.min(), color='blue')
            ax1.plot(dtime[s,:]*2*numpy.pi, data[s,:].imag - data[s,:].imag.min(), color='green')
            
            #ninety = dtime[s, numpy.where( data[s,:].real == data[s,:].real.max() )[0]].mean()
            #ax1.plot(dtime[s,0:200]*2*numpy.pi, (numpy.cos((dtime[s,0:200]-ninety+0.75)*2*numpy.pi)+1)/2*(data[s,:].real.max()-data[s,:].real.min()), color='cyan')
            
            #ninety = dtime[s, numpy.where( data[s,:].imag == data[s,:].imag.max() )[0]].mean()
            #ax1.plot(dtime[s,0:200]*2*numpy.pi, (numpy.cos((dtime[s,0:200]-ninety+0.75)*2*numpy.pi)+1)/2*(data[s,:].imag.max()-data[s,:].imag.min()), color='magenta')
            
            ax1.set_title('Beam %i, Tune. %i, Pol. %i' % (standMapper[s]//4+1, standMapper[s]%4//2+1, standMapper[s]%2))
    else:
        sortedMapper = sorted(standMapper)
        for k, aStand in enumerate(sortedMapper):
            s = standMapper.index(aStand)
            if s+1 not in args.keep:
                continue
            
            ax1 = fig.add_subplot(2, 2, k+1)
            ax1.plot(dtime[s,0:endPt], data[s,0:endPt].real, marker='+', linestyle=' ', color='blue')
            ax1.plot(dtime[s,0:endPt], data[s,0:endPt].imag, marker='+', linestyle=' ', color='green')
            
            ax1.set_title('Beam %i, Tune. %i, Pol. %i' % (standMapper[s]//4+1, standMapper[s]%4//2+1, standMapper[s]%2))
        
    plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="create an eye diagram for some portion of a TBN or DRX file",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
    parser.add_argument('filename', type=str, 
                        help='filename to analyze')
    parser.add_argument('-f', '--freq', type=aph.positive_or_zero_float, default=38.0,
                        help='frequency of the observations in MHz')
    parser.add_argument('-i', '--input-freq', type=aph.positive_or_zero_float, default=38.25,
                        help='frequency of the input sinusoid in MHz')
    parser.add_argument('-s', '--skip', type=aph.positive_or_zero_float, default=0.0,
                        help='number of seconds to skip at the beginning of the file')
    parser.add_argument('-t', '--time', type=aph.positive_float, default=10.0,
                        help='amount of data to plot in seconds')
    parser.add_argument('-k', '--keep', type=aph.csv_int_list, default='1,2,3,4',
                        help='data array indiece (stands/beams) to keep')
    parser.add_argument('-r', '--rectilinear', dest='polar', action='store_false',
                        help='do not plot the eye diagram in polar coordinates')
    args = parser.parse_args()
    args.freq *= 1e6
    args.input_freq *= 1e6
    main(args)
    
    
