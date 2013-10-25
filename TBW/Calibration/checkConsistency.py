#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import numpy

from lsl.imaging import utils
from lsl.astro import DJD_OFFSET
from lsl.statistics import robust

from matplotlib import pyplot as plt


def main(args):
	# Get the list of .cs filenames to parse
	filenames = args
	
	# Go!
	utcs = []
	lsts = []
	delaysX = {}
	delaysY = {}
	order = []
	for filename in filenames:
		## Figure out what the corresponding FITS IDI file is so that we can pull 
		## out the date and LST
		fitsname = filename.replace('.sc', '.FITS_1')
		idi = utils.CorrelatedData(fitsname)
		lo = idi.getObserver()
		lo.date = idi.dateObs.strftime("%Y/%m/%d %H:%M:%S")
		lst = float(lo.sidereal_time()) * 12.0/numpy.pi
		utcs.append(lo.date + DJD_OFFSET)
		lsts.append(lst)
		
		## Load in the actual file
		data = numpy.loadtxt(filename)
		for i in xrange(data.shape[0]):
			stand, ax, dx, ay, dy = data[i,:]
			stand = int(stand)
			dx = float(dx)
			dy = float(dy)
			if stand not in order:
				order.append(stand)
				
			try:
				delaysX[stand].append(dx)
				delaysY[stand].append(dy)
			except KeyError:
				delaysX[stand] = [dx,]
				delaysY[stand] = [dy,]
				
	# Convert to arrays
	for stand in delaysX.keys():
		delaysX[stand] = numpy.array(delaysX[stand])
		delaysY[stand] = numpy.array(delaysY[stand])
		
	# Calculate the mean delay for each capture
	vsX = numpy.zeros((len(delaysX.keys()), len(filenames)))
	vsY = numpy.zeros_like(vsX)
	for i,stand in enumerate(delaysX.keys()):
		dx = delaysX[stand]
		dy = delaysY[stand]
		for j in xrange(len(dx)):
			vsX[i,j] = dx[j]
			vsY[i,j] = dy[j]
	msX = vsX.mean(axis=0)
	msY = vsY.mean(axis=0)
	
	# Merge the delays together
	fh = open('merged.delays', 'w')
	for stand in order:
		delayX = robust.mean( delaysX[stand] - msX)
		dstdX  = robust.std(  delaysX[stand] - msX)
		delayY = robust.mean( delaysY[stand] - msY)
		dstdY  = robust.std(  delaysY[stand] - msY)
		fh.write("%3i  %.2f  %.2f  %.2f  %.2f\n" % (stand, 0.0, delayX, 0.0, delayY))
	fh.close()
	
	#
	# By LST
	#
	lsts = numpy.array(lsts)
	orderL = numpy.argsort(lsts)
	
	figL = plt.figure()
	axLX1 = figL.add_subplot(3, 2, 1)
	axLY1 = figL.add_subplot(3, 2, 2)
	axLX2 = figL.add_subplot(3, 2, 3)
	axLY2 = figL.add_subplot(3, 2, 4)
	axLX3 = figL.add_subplot(3, 2, 5)
	axLY3 = figL.add_subplot(3, 2, 6)
	for stand in delaysX.keys():
		axLX1.plot(lsts[orderL], delaysX[stand][orderL])
		axLY1.plot(lsts[orderL], delaysY[stand][orderL])
	axLX1.set_title('X pol.')
	axLX1.set_xlabel("LST [hours]")
	axLX1.set_ylabel("$\\tau_X$ [ns]")
	axLY1.set_title('Y pol.')
	axLY1.set_xlabel("LST [hours]")
	axLY1.set_ylabel("$\\tau_Y$ [ns]")
	
	axLX2.plot(lsts[orderL], msX[orderL], linestyle='-', marker='x')
	axLY2.plot(lsts[orderL], msY[orderL], linestyle='-', marker='x')
	axLX2.set_xlabel("LST [hours]")
	axLX2.set_ylabel("$|\\tau_x|$ [ns]")
	axLY2.set_xlabel("LST [hours]")
	axLY2.set_ylabel("$|\\tau_y|$ [ns]")
	
	for stand in delaysX.keys():
		dx = delaysX[stand] - msX
		dy = delaysY[stand] - msY
		axLX3.plot(lsts[orderL], dx[orderL])
		axLY3.plot(lsts[orderL], dy[orderL])
	axLX3.set_xlabel("LST [hours]")
	axLX3.set_ylabel("$\\tau_X-|\\tau_x|$ [ns]")
	axLY3.set_xlabel("LST [hours]")
	axLY3.set_ylabel("$\\tau_Y-|\\tau_y|$ [ns]")
	
	figL.tight_layout()
	
	#
	# By JD
	#
	utcs = numpy.array(utcs)
	orderJ = numpy.argsort(utcs)
	utcOffset = utcs.min()
	utcs -= utcOffset
	
	figJ = plt.figure()
	axJX1 = figJ.add_subplot(3, 2, 1)
	axJY1 = figJ.add_subplot(3, 2, 2)
	axJX2 = figJ.add_subplot(3, 2, 3)
	axJY2 = figJ.add_subplot(3, 2, 4)
	axJX3 = figJ.add_subplot(3, 2, 5)
	axJY3 = figJ.add_subplot(3, 2, 6)
	for stand in delaysX.keys():
		axJX1.plot(utcs[orderJ], delaysX[stand][orderJ])
		axJY1.plot(utcs[orderJ], delaysY[stand][orderJ])
	axJX1.set_title('X pol.')
	axJX1.set_xlabel("JD [days - %.1f]" % utcOffset)
	axJX1.set_ylabel("$\\tau_X$ [ns]")
	axJY1.set_title('Y pol.')
	axJY1.set_xlabel("JD [days - %.1f]" % utcOffset)
	axJY1.set_ylabel("$\\tau_Y$ [ns]")
	
	axJX2.plot(utcs[orderJ], msX[orderJ], linestyle='-', marker='x')
	axJY2.plot(utcs[orderJ], msY[orderJ], linestyle='-', marker='x')
	axJX2.set_xlabel("JD [days - %.1f]" % utcOffset)
	axJX2.set_ylabel("$|\\tau_x|$ [ns]")
	axJY2.set_xlabel("JD [days - %.1f]" % utcOffset)
	axJY2.set_ylabel("$|\\tau_y|$ [ns]")
	
	tooHighX = []
	tooHighY = []	
	for stand in delaysX.keys():
		dx = delaysX[stand] - msX
		dy = delaysY[stand] - msY
		print "%3i:  %6.3f +/- %5.3f  %6.3f +/- %5.3f" % (stand, dx.mean(), dx.std(), dy.mean(), dy.std())
		if dx.std() > 1:
			tooHighX.append(stand)
		if dy.std() > 1:
			tooHighY.append(stand)
		axJX3.plot(utcs[orderJ], dx[orderJ])
		axJY3.plot(utcs[orderJ], dy[orderJ])
	axJX3.set_xlabel("JD [days - %.1f]" % utcOffset)
	axJX3.set_ylabel("$\\tau_X-|\\tau_x|$ [ns]")
	axJY3.set_xlabel("JD [days - %.1f]" % utcOffset)
	axJY3.set_ylabel("$\\tau_Y-|\\tau_y|$ [ns]")
	
	print len(tooHighX), len(tooHighY)
	print tooHighX, tooHighY
	
	figJ.tight_layout()
	
	plt.show()


if __name__ == "__main__":
	main(sys.argv[1:])
	
