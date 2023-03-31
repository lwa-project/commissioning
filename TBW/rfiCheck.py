#!/usr/bin/env python3

"""
Script to take a single TBW capture and create a RFI-centered HDF5 file for stands 1, 10, 54, 
248, 251, and 258 (the outlier).  These stands correspond to the four corners of the array, the
center, and the outlier.  The HDF5 contains values for the spectral kurtosis estimated from
the data and various statistics about the timeseries (mean, std. dev., percentiles, etc.)
"""

import os
import sys
import math
import h5py
import numpy
import argparse

from scipy.stats import scoreatpercentile as percentile

from lsl.common import stations
from lsl.reader import tbw, tbn
from lsl.reader import errors
from lsl.correlator import fx as fxc
from lsl.correlator._core import FEngine
from lsl.astro import unix_to_utcjd, DJD_OFFSET
from lsl.common.progress import ProgressBar
from lsl.statistics import kurtosis
from lsl.common.paths import DATA as dataPath
from lsl.misc import parser as aph

import matplotlib.pyplot as plt


def expandMask(mask, radius=2, merge=False):
    """
    Expand a 2-D numpy mask array around existing mask elements.
    """
    
    mask2 = numpy.zeros(mask.shape, dtype=numpy.int16)
    for i in range(mask.shape[0]):
        for j in range(mask.shape[1]):
            if mask[i,j] == 1:
                for k in range(-radius,radius+1):
                    try:
                        mask2[i,j+k] = True
                    except IndexError:
                        pass
                    
    if merge:
        mask3 = mask2.sum(axis=0)
        for i in range(mask.shape[1]):
            if mask3[i] > 0:
                mask2[:,i] = True
    
    mask2 = mask2.astype(numpy.bool)
    
    return mask2


def main(args):
    # Set the station
    if args.metadata is not None:
        station = stations.parse_ssmif(args.metadata)
        ssmifContents = open(args.metadata).readlines()
    else:
        station = stations.lwa1
        ssmifContents = open(os.path.join(dataPath, 'lwa1-ssmif.txt')).readlines()
    antennas = station.antennas
    
    toKeep = []
    for g in (1, 10, 54, 248, 251, 258):
        for i,ant in enumerate(antennas):
            if ant.stand.id == g and ant.pol == 0:
                toKeep.append(i)
    for i,j in enumerate(toKeep):
        print(i, j, antennas[j].stand.id)

    # Length of the FFT
    LFFT = args.fft_length

    # Make sure that the file chunk size contains is an integer multiple
    # of the FFT length so that no data gets dropped
    maxFrames = int((30000*260)/float(LFFT))*LFFT
    # It seems like that would be a good idea, however...  TBW data comes one
    # capture at a time so doing something like this actually truncates data 
    # from the last set of stands for the first integration.  So, we really 
    # should stick with
    maxFrames = (30000*260)

    fh = open(args.filename, "rb")
    nFrames = os.path.getsize(args.filename) // tbw.FRAME_SIZE
    dataBits = tbw.get_data_bits(fh)
    # The number of ant/pols in the file is hard coded because I cannot figure out 
    # a way to get this number in a systematic fashion
    antpols = len(antennas)
    nChunks = int(math.ceil(1.0*nFrames/maxFrames))
    if dataBits == 12:
        nSamples = 400
    else:
        nSamples = 1200

    # Read in the first frame and get the date/time of the first sample 
    # of the frame.  This is needed to get the list of stands.
    junkFrame = tbw.read_frame(fh)
    fh.seek(0)
    beginTime = junkFrame.time
    beginDate = junkFrame.time.datetime

    # File summary
    print("Filename: %s" % args.filename)
    print("Date of First Frame: %s" % str(beginDate))
    print("Ant/Pols: %i" % antpols)
    print("Sample Length: %i-bit" % dataBits)
    print("Frames: %i" % nFrames)
    print("Chunks: %i" % nChunks)
    print("===")

    nChunks = 1

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

    # Master loop over all of the file chunks
    masterSpectra = numpy.zeros((nChunks, antpols, LFFT))
    for i in range(nChunks):
        # Find out how many frames remain in the file.  If this number is larger
        # than the maximum of frames we can work with at a time (maxFrames),
        # only deal with that chunk
        framesRemaining = nFrames - i*maxFrames
        if framesRemaining > maxFrames:
            framesWork = maxFrames
        else:
            framesWork = framesRemaining
        print("Working on chunk %i, %i frames remaining" % ((i+1), framesRemaining))

        data = numpy.zeros((12, 12000000), dtype=numpy.int16)
        # If there are fewer frames than we need to fill an FFT, skip this chunk
        if data.shape[1] < 2*LFFT:
            break
        # Inner loop that actually reads the frames into the data array
        for j in range(framesWork):
            # Read in the next frame and anticipate any problems that could occur
            try:
                cFrame = tbw.read_frame(fh)
            except errors.EOFError:
                break
            except errors.SyncError:
                #print("WARNING: Mark 5C sync error on frame #%i" % (int(fh.tell())/tbw.FRAME_SIZE-1))
                continue
            if not cFrame.header.is_tbw:
                continue
            
            stand = cFrame.header.id
            # In the current configuration, stands start at 1 and go up to 10.  So, we
            # can use this little trick to populate the data array
            aStand = 2*(stand-1)
            #if cFrame.header.frame_count % 10000 == 0 and config['verbose']:
                #print("%3i -> %3i  %6.3f  %5i  %i" % (stand, aStand, cFrame.time, cFrame.header.frame_count, cFrame.payload.timetag))

            # Actually load the data.  x pol goes into the even numbers, y pol into the 
            # odd numbers
            count = cFrame.header.frame_count - 1
            if aStand not in toKeep:
                continue
            
            # Convert to reduced index
            aStand = 2*toKeep.index(aStand)
            
            data[aStand,   count*nSamples:(count+1)*nSamples] = cFrame.payload.data[0,:]
            data[aStand+1, count*nSamples:(count+1)*nSamples] = cFrame.payload.data[1,:]
    
        # Time series analysis - mean, std. dev, saturation count
        tsMean = data.mean(axis=1)
        tsStd = data.std(axis=1)
        tsSat = numpy.where( (data == 2047) | (data == -2047), 1, 0 ).sum(axis=1)
        
        # Time series analysis - percentiles
        p = [50, 75, 90, 95, 99]
        tsPct = numpy.zeros((data.shape[0], len(p)))
        for i in range(len(p)):
            for j in range(data.shape[0]):
                tsPct[j,i] = percentile(numpy.abs(data[j,:]), p[i])
    
        # Frequency domain analysis - spectra
        freq = numpy.fft.fftfreq(2*args.fft_length, d=1.0/196e6)
        freq = freq[:args.fft_length]
        
        delays = numpy.zeros((data.shape[0], freq.size))
        signalsF, validF = FEngine(data, freq, delays, LFFT=args.fft_length, Overlap=1, sample_rate=196e6, clip_level=0)
        
        # Cleanup to save memory
        del validF, data
        print(signalsF.shape)
        
        # SK control values
        skM = signalsF.shape[2]
        skN = 1
        
        # Frequency domain analysis -  spectral kurtosis
        k = numpy.zeros((signalsF.shape[0], signalsF.shape[1]))
        for l in range(signalsF.shape[0]):
            for m in range(freq.size):
                k[l,m] = kurtosis.spectral_fft(signalsF[l,m,:])
        kl, kh = kurtosis.get_limits(4, skM, skN)
        print(kl, kh)
        
        # Integrate the spectra for as long as we can
        masterSpectra = (numpy.abs(signalsF)**2).mean(axis=2)
        del signalsF
        
        # Mask out bad values (high spectral kurtosis) for the plot
        mask = numpy.where( (k < kl) | ( k > kh), 1, 0 )
        mask = expandMask(mask, radius=4, merge=True)
        
        masterSpectra = numpy.ma.array(masterSpectra, mask=mask)
        
        # Save the data to an HDF5 file
        outname = os.path.splitext(args.filename)[0]
        outname = "%s-RFI.hdf5" % outname
        
        f = h5py.File(outname, 'w')
        f.attrs['filename'] = args.filename
        f.attrs['mode'] = 'TBW'
        f.attrs['station'] = 'LWA-1'
        f.attrs['dataBits'] = dataBits
        f.attrs['startTime'] = beginTime
        f.attrs['startTime_units'] = 's'
        f.attrs['startTime_sys'] = 'unix'
        f.attrs['sample_rate'] = 196e6
        f.attrs['sample_rate_units'] = 'Hz'
        f.attrs['RBW'] = freq[1]-freq[0]
        f.attrs['RBW_Units'] = 'Hz'
        
        f.attrs['SK-M'] = skM
        f.attrs['SK-N'] = skN
        
        for l in range(len(toKeep)):
            antX = antennas[toKeep[l]]
            antY = antennas[toKeep[l]+1]
            
            stand = f.create_group('Stand%03i' % antX.stand.id)
            stand['freq'] = freq
            stand['freq'].attrs['Units'] = 'Hz'
            
            polX = stand.create_group('X')
            polY = stand.create_group('Y')
            polX.attrs['tsMean'] = tsMean[2*l]
            polY.attrs['tsMean'] = tsMean[2*l+1]
            polX.attrs['tsStd'] = tsStd[2*l]
            polY.attrs['tsStd'] = tsStd[2*l+1]
            polX.attrs['tsSat'] = tsSat[2*l]
            polY.attrs['tsSat'] = tsSat[2*l+1]
            for i,v in enumerate(p):
                polX.attrs['ts%02i' % v] = tsPct[2*l][i]
                polY.attrs['ts%02i' % v] = tsPct[2*l+1][i]
            
            polX['spectrum'] = masterSpectra[2*l,:]
            polX['spectrum'].attrs['axis0'] = 'frequency'
            polY['spectrum'] = masterSpectra[2*l+1,:]
            polY['spectrum'].attrs['axis0'] = 'frequency'
            
            polX['kurtosis'] = k[2*l,:]
            polX['kurtosis'].attrs['axis0'] = 'frequency'
            polY['kurtosis'] = k[2*l+1,:]
            polY['kurtosis'].attrs['axis0'] = 'frequency'
        
        # The plot
        fig = plt.figure()
        ax1 = fig.add_subplot(2, 1, 1)
        ax2 = fig.add_subplot(2, 1, 2)
        for l in range(k.shape[0]):
            ant = antennas[toKeep[l/2]]
            
            ax1.plot(freq/1e6, numpy.log10(masterSpectra[l,:])*10, label='Stand %i, Pol %i' % (ant.stand.id, ant.pol+l%2))
            
            ax2.plot(freq/1e6, k[l,:], label='Stand %i, Pol %i' % (ant.stand.id, ant.pol+l%2))
            
        ax2.hlines(kl, freq[0]/1e6, freq[-1]/1e6, linestyle=':', label='Kurtosis Limit 4$\sigma$')
        ax2.hlines(kh, freq[0]/1e6, freq[-1]/1e6, linestyle=':', label='Kurtosis Limit 4$\sigma$')
        
        ax1.set_xlabel('Frequency [MHz]')
        ax1.set_ylabel('PSD [arb. dB/RBW]')
        ax1.legend(loc=0)
        
        ax2.set_ylim((kl/2, kh*2))
        ax2.set_xlabel('Frequency [MHz]')
        ax2.set_ylabel('Spectral Kurtosis')
        ax2.legend(loc=0)
        
        plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='read in TBW files and create a collection of RFI statistics', 
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
    parser.add_argument('filename', type=str, nargs='+', 
                        help='filename to process')
    parser.add_argument('-m', '--metadata', type=str, 
                        help='name of the SSMIF or metadata tarball file to use for mappings')
    parser.add_argument('-l', '--fft-length', type=aph.positive_int, default=4096, 
                        help='set FFT length')
    args = parser.parse_args()
    main(args)
    
