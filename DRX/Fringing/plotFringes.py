#!/usr/bin/env python

"""
Simple script to plot up the NPZ files created by fringeDipole.py/fringeBeam.py.
"""

# Python3 compatiability
from __future__ import print_function, division
import sys
if sys.version_info > (3,):
    xrange = range
    
import os
import sys
import glob
import numpy
import argparse

from datetime import datetime

from matplotlib import pyplot as plt


def main(args):
    times = []
    amp1 = []
    amp2 = []
    amp3 = []
    amp4 = []
    amp5 = []
    amp6 = []
    phs1 = []
    phs2 = []

    for filename in args.filenames:
        dataDict = numpy.load(filename)

        srate = dataDict['srate']
        tStart = dataDict['tStart']
        tStart = datetime.utcfromtimestamp(tStart)
        
        freq1 = dataDict['freq1']
        vis1 = dataDict['vis1'][1,freq1.size//4:freq1.size*3//4]
        auto11 = dataDict['vis1'][0,freq1.size//4:freq1.size*3//4]
        auto12 = dataDict['vis1'][2,freq1.size//4:freq1.size*3//4]
        
        freq2 = dataDict['freq2']
        vis2 = dataDict['vis2'][1,freq2.size//4:freq2.size*3//4]
        auto21 = dataDict['vis2'][0,freq2.size//4:freq2.size*3//4]
        auto22 = dataDict['vis2'][2,freq2.size//4:freq2.size*3//4]

        times.append( tStart)
        amp1.append( numpy.abs(vis1).mean() )
        amp2.append( numpy.abs(vis2).mean() )
        amp3.append( numpy.abs(auto11).mean() )
        amp4.append( numpy.abs(auto21).mean() )
        amp5.append( numpy.abs(auto12).mean() )
        amp6.append( numpy.abs(auto22).mean() )
        phs1.append( numpy.angle(vis1).mean() )
        phs2.append( numpy.angle(vis2).mean() )
        
        dataDict.close()

    amp1 = numpy.array(amp1)
    amp2 = numpy.array(amp2)
    amp3 = numpy.array(amp3)
    amp4 = numpy.array(amp4)
    amp5 = numpy.array(amp5)
    amp6 = numpy.array(amp6)
    phs1 = numpy.array(phs1)
    phs2 = numpy.array(phs2)

    #
    # Amplitude
    #

    fig = plt.figure()
    ax1 = fig.add_subplot(1, 2, 1)
    ax2 = fig.add_subplot(1, 2, 2)
    fig.autofmt_xdate()
    ax1b = ax1.twinx()
    ax2b = ax2.twinx()

    l1, = ax1.plot_date(times, amp1, linestyle='-', label='Beam-Dipole')
    l2, = ax1b.plot_date(times, amp3, linestyle='--', color='green', alpha=0.40, label='Beam-Beam')
    l3, = ax2.plot_date(times, amp2, linestyle='-', label='Beam-Dipole')
    l4, = ax2b.plot_date(times, amp4, linestyle='--', color='green', alpha=0.40, label='Beam-Beam')

    fig.suptitle("%s to %s UTC" % (times[0].strftime("%Y/%m/%d %H:%M"), times[-1].strftime("%Y/%m/%d %H:%M")))
    ax1.set_xlabel('Time')
    ax2.set_xlabel('Time')
    ax1.set_ylabel('Vis. Amp. [arb.]')
    ax2.set_ylabel('Vis. Amp. [arb.]')
    ax1.set_title('%.1f MHz @ %.2f MHz BW' % (freq1.mean()/1e6, 0.75*srate/1e6))
    ax2.set_title('%.1f MHz @ %.2f MHz BW' % (freq2.mean()/1e6, 0.75*srate/1e6))
    ax1.legend([l1, l2], ['Beam-Outlier', 'Beam-Beam'], loc=0)
    ax2.legend([l3, l4], ['Beam-Outlier', 'Beam-Beam'], loc=0)

    r1 = amp1.max() - amp2.min()
    r2 = amp2.max() - amp2.min()
    r3 = amp3.max() - amp3.min()
    r4 = amp4.max() - amp4.min()
    ax1.set_ylim((amp1.min() - 0.05*r1, amp1.max() + 0.05*r1))
    ax2.set_ylim((amp2.min() - 0.05*r2, amp2.max() + 0.05*r2))
    ax1b.set_ylim((amp3.min() - 0.05*r3, amp3.max() + 0.05*r3))
    ax2b.set_ylim((amp4.min() - 0.05*r4, amp4.max() + 0.05*r4))
    
    plt.draw()
    
    #
    # Phase
    #
    
    fig = plt.figure()
    ax1 = fig.add_subplot(1, 2, 1)
    ax2 = fig.add_subplot(1, 2, 2)

    ax1.plot_date(times, phs1*180/numpy.pi, linestyle='-')
    ax2.plot_date(times, phs2*180/numpy.pi, linestyle='-')

    fig.suptitle("%s to %s UTC" % (times[0].strftime("%Y/%m/%d %H:%M"), times[-1].strftime("%Y/%m/%d %H:%M")))
    ax1.set_xlabel('Time')
    ax2.set_xlabel('Time')
    ax1.set_ylabel('Vis. Phase [deg.]')
    ax2.set_ylabel('Vis. Phase [deg.]')
    ax1.set_title('%.1f MHz @ %.2f MHz BW' % (freq1.mean()/1e6, 0.75*srate/1e6))
    ax2.set_title('%.1f MHz @ %.2f MHz BW' % (freq2.mean()/1e6, 0.75*srate/1e6))

    fig.autofmt_xdate()
    
    plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Simple script to plot up the NPZ files created by fringeDipole.py/fringeBeam.py',
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('filenames', nargs='+',
            help='NPZ files to plot.')

    args = parser.parse_args()
    main(args)

