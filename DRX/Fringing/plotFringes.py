#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Simple script to plot up the NPZ files created by fringeDipole.py/fringeBeam.py.
"""

import os
import sys
import glob
import numpy

from datetime import datetime

from matplotlib import pyplot as plt


def main(args):
	fig = plt.figure()
	ax1 = fig.add_subplot(1, 2, 1)
	ax2 = fig.add_subplot(1, 2, 2)

	times = []
	amp1 = []
	amp2 = []
	amp3 = []
	amp4 = []

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

	amp1 = numpy.array(amp1)
	amp2 = numpy.array(amp2)
	amp3 = numpy.array(amp3)
	amp4 = numpy.array(amp4)

	ax1.plot_date(times, amp1, linestyle='-')
	#ax1.plot_date(times, amp3*amp1.mean()/amp3.mean(), linestyle='--')
	ax2.plot_date(times, amp2, linestyle='-')
	#ax2.plot_date(times, amp4*amp2.mean()/amp4.mean(), linestyle='--')

	ax1.set_xlabel('Time')
	ax2.set_xlabel('Time')
	ax1.set_ylabel('Vis. Amp. [arb.]')
	ax2.set_ylabel('Vis. Amp. [arb.]')
	ax1.set_title('%.1f MHz @ %.2f MHz BW' % (freq1.mean()/1e6, 0.75*srate/1e6))
	ax2.set_title('%.1f MHz @ %.2f MHz BW' % (freq2.mean()/1e6, 0.75*srate/1e6))

	fig.autofmt_xdate()
	plt.show()


if __name__ == "__main__":
	main(sys.argv[1:])

