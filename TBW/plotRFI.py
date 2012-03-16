#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Plot the output of rfiCheck.py.

$Rev$
$LastChangedBy$
$LastChangedDate$
"""

import os
import sys
import h5py
import numpy
from datetime import datetime

from lsl.misc.mathutil import to_dB
from lsl.statistics import kurtosis

from matplotlib import pyplot as plt


def main(args):
	filename = args[0]
	h = h5py.File(filename, 'r')
	tStart = datetime.utcfromtimestamp(h.attrs['startTime'])
	
	print "Filename: %s" % filename
	print "Date: %s" % tStart
	
	print "Groups Found:"
	stands = list(h)
	for s in stands:
		sn = int(s[-3:])
		print "  #%i" % sn
	print " "
	
	for s in stands:
		sn = int(s[-3:])
		print "Stand #%i" % sn
		
		stand = h.get(s, None)
		polX = stand.get('X', None)
		polY = stand.get('Y', None)
		
		# Get timeseries statistics
		tMeanX = polX.attrs['tsMean']
		tMeanY = polY.attrs['tsMean']
		tStdX  = polX.attrs['tsStd']
		tStdY  = polY.attrs['tsStd']
		tSatX  = polX.attrs['tsSat']
		tSatY  = polY.attrs['tsSat']
		t90X   = polX.attrs['ts90']
		t90Y   = polY.attrs['ts90']
		
		print "  X Pol. Timeseries Statistics"
		print "    Mean:          %.3f" % tMeanX
		print "    Std. Dev.:     %.3f" % tStdX
		print "    Saturation:    %i" % tSatX
		print "    90-percentile: %i" % t90X 
		print "  Y Pol. Timeseries Statistics"
		print "    Mean:          %.3f" % tMeanY
		print "    Std. Dev.:     %.3f" % tStdY
		print "    Saturation:    %i" % tSatY
		print "    90-percentile: %i" % t90Y
		
		# Get frequency
		freq = numpy.zeros(stand['freq'].shape, dtype=stand['freq'].dtype)
		stand['freq'].read_direct(freq)
		
		# Spectra - X & Y
		specX = numpy.zeros(polX['spectrum'].shape, dtype=polX['spectrum'].dtype)
		polX['spectrum'].read_direct(specX)
		specY = numpy.zeros(polY['spectrum'].shape, dtype=polY['spectrum'].dtype)
		polY['spectrum'].read_direct(specY)
		
		# Kurtosis - X & Y
		skX = numpy.zeros(polX['kurtosis'].shape, dtype=polX['kurtosis'].dtype)
		polX['kurtosis'].read_direct(skX)
		skY = numpy.zeros(polY['kurtosis'].shape, dtype=polY['kurtosis'].dtype)
		polY['kurtosis'].read_direct(skY)
		
		# 4sigma kurtosis limits
		kl, kh = kurtosis.getLimits(4, h.attrs['SK-M'], h.attrs['SK-N'])
		
		fig = plt.figure()
		ax1 = fig.add_subplot(2, 1, 1)
		ax2 = fig.add_subplot(2, 1, 2)
		ax1.plot(freq/1e6, to_dB(specX), label='X')
		ax1.plot(freq/1e6, to_dB(specY), label='Y')
		
		ax1.set_xlim((5,95))
		ax1.set_xlabel('Frequency [MHz]')
		ax1.set_ylabel('Power [arb. dB/RBW]')
		ax1.legend(loc=0)
		
		ax2.plot(freq/1e6, skX, label='X')
		ax2.plot(freq/1e6, skY, label='Y')
		ax2.hlines(kl, freq[0]/1e6, freq[-1]/1e6, linestyle=':', label='4$\sigma$ Limits')
		ax2.hlines(kh, freq[0]/1e6, freq[-1]/1e6, linestyle=':')
		
		ax2.set_xlim((5,95))
		ax2.set_ylim((kl/4,kh*4))
		ax2.set_xlabel('Frequency [MHz]')
		ax2.set_ylabel('Spectral Kurtosis')
		handles, labels = ax2.get_legend_handles_labels()
		ax2.legend(handles[:-1], labels[:-1], loc=0)
		
		fig.suptitle('Stand #%i' % sn)
		print " "
	
	plt.show()


if __name__ == "__main__":
	main(sys.argv[1:])
	
