#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Simple script to plot up the NPZ files created by fringeDipole.py/fringeBeam.py.

$Rev$
$LastChangedBy$
$LastChangedDate$
"""

import os
import sys
import glob
import numpy

from datetime import datetime

from matplotlib import pyplot as plt


def main(args):
	times = []
	amp1 = []
	amp2 = []
	amp3 = []
	amp4 = []
	phs1 = []
	phs2 = []

	for filename in args:
		dataDict = numpy.load(filename)

		srate = dataDict['srate']
		tStart = datetime.utcfromtimestamp(dataDict['tStart'])
		
		freq1 = dataDict['freq1']
		vis1 = dataDict['vis1'][1,freq1.size/4:freq1.size*3/4]
		auto1 = dataDict['vis1'][0,freq1.size/4:freq1.size*3/4]
		
		freq2 = dataDict['freq2']
		vis2 = dataDict['vis2'][1,freq1.size/4:freq1.size*3/4]
		auto2 = dataDict['vis2'][0,freq1.size/4:freq1.size*3/4]

		times.append( tStart)
		amp1.append( numpy.abs(vis1).mean() )
		amp2.append( numpy.abs(vis2).mean() )
		amp3.append( numpy.abs(auto1).mean() )
		amp4.append( numpy.abs(auto2).mean() )
		phs1.append( numpy.angle(vis1).mean() )
		phs2.append( numpy.angle(vis2).mean() )

	amp1 = numpy.array(amp1)
	amp2 = numpy.array(amp2)
	amp3 = numpy.array(amp3)
	amp4 = numpy.array(amp4)
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

	ax1.plot_date(times, amp1, linestyle='-')
	ax1b.plot_date(times, amp3, linestyle='--', color='green', alpha=0.40)
	ax2.plot_date(times, amp2, linestyle='-')
	ax2b.plot_date(times, amp4, linestyle='--', color='green', alpha=0.40)

	fig.suptitle("%s to %s UTC" % (times[0].strftime("%Y/%m/%d %H:%M"), times[-1].strftime("%Y/%m/%d %H:%M")))
	ax1.set_xlabel('Time')
	ax2.set_xlabel('Time')
	ax1.set_ylabel('Vis. Amp. [arb.]')
	ax2.set_ylabel('Vis. Amp. [arb.]')
	ax1.set_title('%.1f MHz @ %.2f MHz BW' % (freq1.mean()/1e6, 0.75*srate/1e6))
	ax2.set_title('%.1f MHz @ %.2f MHz BW' % (freq2.mean()/1e6, 0.75*srate/1e6))

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

