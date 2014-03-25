#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Simple script to plot up the NPZ files created by fringeDipole.py/fringeBeam.py.

$Rev: 1244 $
$LastChangedBy: jayce $
$LastChangedDate: 2013-03-14 10:43:10 -0600 (Thu, 14 Mar 2013) $
"""

import os
import sys
import glob
import numpy
import h5py

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
    
    infile = h5py.File(args[0],'r')
    srate = infile.attrs['SRATE']
    
    # Get timestamps
    time = numpy.zeros(infile['Time']['Timesteps'].shape, dtype=infile['Time']['Timesteps'].dtype)
    infile['Time']['Timesteps'].read_direct(time)
    
    # Get Freqs
    freq1 = numpy.zeros(infile['Frequencies']['Tuning1'].shape, dtype=infile['Frequencies']['Tuning1'].dtype)
    freq2 = numpy.zeros(infile['Frequencies']['Tuning2'].shape, dtype=infile['Frequencies']['Tuning2'].dtype)
    infile['Frequencies']['Tuning1'].read_direct(freq1)
    infile['Frequencies']['Tuning2'].read_direct(freq2)
    
    for i in range(time.shape[0]):
        times.append(datetime.utcfromtimestamp(time[i][0]))
        
        vis1 = infile['Visibilities']['Tuning1'][i, 1,freq1.size/4:freq1.size*3/4]
        auto11 = infile['Visibilities']['Tuning1'][i, 0,freq1.size/4:freq1.size*3/4]
        auto12 = infile['Visibilities']['Tuning1'][i, 2,freq1.size/4:freq1.size*3/4]
        
        vis2 = infile['Visibilities']['Tuning2'][i, 1,freq2.size/4:freq2.size*3/4]
        auto21 = infile['Visibilities']['Tuning2'][i, 0,freq2.size/4:freq2.size*3/4]
        auto22 = infile['Visibilities']['Tuning2'][i, 2,freq2.size/4:freq2.size*3/4]
    
        amp1.append( numpy.abs(vis1).mean() )
        amp2.append( numpy.abs(vis2).mean() )
        amp3.append( numpy.abs(auto11).mean() )
        amp4.append( numpy.abs(auto21).mean() )
        amp5.append( numpy.abs(auto12).mean() )
        amp6.append( numpy.abs(auto22).mean() )
        phs1.append( numpy.angle(vis1).mean() )
        phs2.append( numpy.angle(vis2).mean() )
        
        del vis1, auto11, auto12,  vis2,  auto21, auto22
    
    amp1 = numpy.array(amp1)
    amp2 = numpy.array(amp2)
    amp3 = numpy.array(amp3)
    amp4 = numpy.array(amp4)
    amp5 = numpy.array(amp5)
    amp6 = numpy.array(amp6)
    phs1 = numpy.array(phs1)
    phs2 = numpy.array(phs2)
    
    infile.close()
    
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
	main(sys.argv[1:])

