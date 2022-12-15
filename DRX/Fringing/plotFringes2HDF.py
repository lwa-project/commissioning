#!/usr/bin/env python3

"""
A fancier version of plotFringesHDF.py that makes waterfall-like plots.
"""

# Python2 compatibility
from __future__ import print_function, division
try:
    range = xrange
except NameError:
    pass
    
import os
import sys
import glob
import numpy
import h5py
import argparse

from datetime import datetime

from lsl.statistics import robust
from lsl.misc.mathutils import to_dB

from matplotlib import pyplot as plt


def spectralKurtosis(x, N=1):
    """
    Compute the spectral kurtosis for a set of power measurments averaged
    over N FFTs.  For a distribution consistent with Gaussian noise, this
    value should be ~1.
    """
    
    M = len(x)
    
    k = M*(x**2).sum()/(x.sum())**2 - 1.0
    k *= (M*N+1)/(M-1)
    
    return k


def main(args):
    filenames = args.filenames
    filenames.sort()
    
    fig1 = plt.figure()
    ax1 = fig1.add_subplot(2, 1, 1)
    ax2 = fig1.add_subplot(2, 1, 2)
    
    fig2 = plt.figure()
    ax3 = fig2.add_subplot(2, 1, 1)
    ax4 = fig2.add_subplot(2, 1, 2)
    
    times = []
    vis1 = []
    vis2 = []
    
    infile = h5py.File(args.filenames[0],'r')
    srate = infile.attrs['SRATE']
    
    # Get timestamps
    time = numpy.zeros(infile['Time']['Timesteps'].shape, dtype=infile['Time']['Timesteps'].dtype)
    infile['Time']['Timesteps'].read_direct(time)
    
    # Get Freqs
    freq1 = numpy.zeros(infile['Frequencies']['Tuning1'].shape, dtype=infile['Frequencies']['Tuning1'].dtype)
    freq2 = numpy.zeros(infile['Frequencies']['Tuning2'].shape, dtype=infile['Frequencies']['Tuning2'].dtype)
    infile['Frequencies']['Tuning1'].read_direct(freq1)
    infile['Frequencies']['Tuning2'].read_direct(freq2)
    
    tInt = time[0][1]
    
    for i in range(time.shape[0]):
        times.append(datetime.utcfromtimestamp(time[i][0]))
        vis1.append(infile['Visibilities']['Tuning1'][i, 1,:])
        vis2.append(infile['Visibilities']['Tuning2'][i, 1,:])
    
    infile.close()

    N = srate*tInt/(len(freq1)+1)
    
    print("Got %i files from %s to %s (%s)" % (len(filenames), times[0].strftime("%Y/%m/%d %H:%M:%S"), times[-1].strftime("%Y/%m/%d %H:%M:%S"), (times[-1]-times[0])))
    
    iTimes = []
    for i in range(1, len(times)):
        dt = times[i] - times[i-1]
        iTimes.append(dt.days*24*3600 + dt.seconds + dt.microseconds/1e6)
    iTimes = numpy.array(iTimes)
    print(" -> Interval: %.3f +/- %.3f seconds (%.3f to %.3f seconds)" % (iTimes.mean(), iTimes.std(), iTimes.min(), iTimes.max()))
    
    print("Number of frequency channels: %i (~%.1f Hz/channel)" % (len(freq1)+1, freq1[1]-freq1[0]))
    
    dTimes = []
    for t in times:
        dTimes.append( (t-times[0]).seconds )
    
    freq1 /= 1e6
    freq2 /= 1e6
    
    vis1 = numpy.array(vis1)
    vis1 = numpy.ma.array(vis1, mask=~numpy.isfinite(vis1))
    vis2 = numpy.array(vis2)
    vis2 = numpy.ma.array(vis2, mask=~numpy.isfinite(vis2))
    
    sk = numpy.zeros(freq1.size)
    for i in range(vis1.shape[1]):
        sk[i] = spectralKurtosis(numpy.abs(vis1[:,i])**2, N=N)
    
    skM = robust.mean(sk)
    skS = robust.std(sk)
    bad = numpy.where( numpy.abs(sk - skM) > 4*skS )[0]
    #vis1.mask[:,bad] = True
    
    sk = numpy.zeros_like(freq2)
    for i in range(vis2.shape[1]):
        sk[i] = spectralKurtosis(numpy.abs(vis2[:,i])**2, N=N)
    
    skM = robust.mean(sk)
    skS = robust.std(sk)
    bad = numpy.where( numpy.abs(sk - skM) > 4*skS )[0]
    #vis2.mask[:,bad] = True
    
    data = 1.0*numpy.abs(vis1)
    data = data.ravel()
    data.sort()
    vmin1 = data[int(round(0.15*len(data)))]
    vmax1 = data[int(round(0.85*len(data)))]
    print('Plot range for tuning 1:', vmin1, vmax1)
    
    data = 1.0*numpy.abs(vis2)
    data = data.ravel()
    data.sort()
    vmin2 = data[int(round(0.15*len(data)))]
    vmax2 = data[int(round(0.85*len(data)))]
    print('Plot range for tuning 2:', vmin2, vmax2)
    
    ax1.imshow(numpy.abs(vis1), extent=(freq1[0], freq1[-1], dTimes[0], dTimes[-1]), origin='lower', vmin=vmin1, vmax=vmax1)
    ax2.imshow(numpy.abs(vis2), extent=(freq2[0], freq2[-1], dTimes[0], dTimes[-1]), origin='lower', vmin=vmin2, vmax=vmax2)
    
    ax1.axis('auto')
    ax2.axis('auto')
    
    fig1.suptitle("%s to %s UTC" % (times[0].strftime("%Y/%m/%d %H:%M"), times[-1].strftime("%Y/%m/%d %H:%M")))
    ax1.set_xlabel('Frequency [MHz]')
    ax2.set_xlabel('Frequency [MHz]')
    ax1.set_ylabel('Elapsed Time [s]')
    ax2.set_ylabel('Elapsed Time [s]')
    
    ax3.plot(freq1, numpy.abs(vis1).mean(axis=0))
    ax4.plot(freq2, numpy.abs(vis2).mean(axis=0))
    
    fig2.suptitle("%s to %s UTC" % (times[0].strftime("%Y/%m/%d %H:%M"), times[-1].strftime("%Y/%m/%d %H:%M")))
    ax3.set_xlabel('Frequency [MHz]')
    ax4.set_xlabel('Frequency [MHz]')
    ax3.set_ylabel('Mean Vis. Amp. [arb.]')
    ax4.set_ylabel('Mean Vis. Amp. [arb.]')
    
    plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='A fancier version of plotFringesHDF.py that makes waterfall-like plots.',
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('filenames', nargs='+',
            help='HDF5 file to plot')

    args = parser.parse_args()
    main(args)
