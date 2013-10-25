# -*- coding: utf-8 -*-

"""Simple self-calibration module for correlated TBW and TBN data.

.. versionadded:: 0.5.5
"""

import sys
import math
import numpy

import ols

from lsl.imaging import utils
from lsl.sim.vis import scaleData

__version__ = '0.1'
__revision__ = '$Rev: 497 $'
__all__ = ['selfCal', '__version__', '__revision__', '__all__']

N_CHAN = 500


def unwrap(x, p):
	total = 0
	s = 1e6
	while abs(s) > 0.01:
		slopes = []
		for i in xrange(x.size):
			for j in xrange(x.size):
				if i == j:
					continue
				elif abs(i-j) > x.size/20:
					continue
				s = (p[j] - p[i]) / (x[j] - x[i])
				slopes.append(s)
		slopes = numpy.array(slopes)

		s = numpy.median(slopes)
		ep = s*x
		ep %= (2*numpy.pi)
		#pylab.plot(x, ep)
		
		p -= ep
		total += s
		print s, total
	
	ep = total*x
	c = numpy.round((p - ep) / (2*numpy.pi))
	return ep


def _buildAmplitudeA(aa, dataDict, simDict, chan, pol, refAnt=0):
	"""
	Build the matrix A.
	"""
	
	# Get the baseline and stand counts
	nBLs = len(dataDict['bls'][pol])
	nStands = len(aa.ants)
	
	A = numpy.zeros((nBLs, nStands))
	for i,(l,m) in enumerate(dataDict['bls'][pol]):
		A[i,l] = 1
		A[i,m] = 1
		
	A = numpy.matrix(A)
	return A


def _buildAmplitudeC(aa, dataDict, simDict, chan, pol, refAnt=0, plot=False):
	# Get the baseline and stand counts
	nBLs = len(dataDict['bls'][pol])
	nStands = len(aa.ants)
	
	# Frequency in GHz so that the delays can be in ns
	fq = dataDict['freq'][chan] / 1e9
	
	obsVis = []
	for vis in dataDict['vis'][pol]:
		obsVis.append( numpy.array(vis[chan]) )
	obsVis = numpy.array(obsVis)
	
	simVis = []
	for vis in simDict['vis'][pol]:
		simVis.append( numpy.array(vis[chan]) )
	simVis = numpy.array(simVis)
	
	C = numpy.log(numpy.abs(simVis)) - numpy.log(numpy.abs(obsVis))
	C = C[:,len(chan)/2]
	
	return C


def _buildPhaseA(aa, dataDict, simDict, chan, pol, refAnt=0):
	"""
	Build the matrix A.
	"""
	
	# Get the baseline and stand counts
	nBLs = len(dataDict['bls'][pol])
	nStands = len(aa.ants)
	
	# Frequency in GHz so that the delays can be in ns
	fq = dataDict['freq'][chan] / 1e9
	fIndexList = range(0, fq.size-1, max([fq.size/N_CHAN,1]))
	if len(fIndexList) < 3:
		fIndexList = [0, fq.size/2, -1]
		
	A = numpy.zeros((len(fIndexList)*nBLs, nStands-1))
	for i,f in enumerate(fIndexList):
		for j,(l,m) in enumerate(dataDict['bls'][pol]):
			if l < refAnt:
				A[j+i*nBLs,l]   =  1
			elif l > refAnt:
				A[j+i*nBLs,l-1] =  1
			else:
				pass
				
			if m < refAnt:
				A[j+i*nBLs,m]   = -1
			elif m > refAnt:
				A[j+i*nBLs,m-1] = -1
			else:
				pass
		
	A = numpy.matrix(A)
	return A


def _buildPhaseC(aa, dataDict, simDict, chan, pol, refAnt=0, plot=False):
	# Get the baseline and stand counts
	nBLs = len(dataDict['bls'][pol])
	nStands = len(aa.ants)
	
	# Frequency in GHz so that the delays can be in ns
	fq = dataDict['freq'][chan] / 1e9
	
	obsVis = []
	for vis in dataDict['vis'][pol]:
		obsVis.append( numpy.array(vis[chan]) )
	obsVis = numpy.array(obsVis)
	
	simVis = []
	for vis in simDict['vis'][pol]:
		simVis.append( numpy.array(vis[chan]) )
	simVis = numpy.array(simVis)
	
	#from test2 import unwrapFFT1D
	C = numpy.angle(simVis / obsVis)
	#for i in xrange(C.shape[0]):
	#	C[i,:] = unwrapFFT1D(fq, simVis[i,:]/obsVis[i,:])
	
	if plot:
		resid = (C[:,[0,fq.size/2,-1]]**2).sum(axis=1)
		o = sorted(resid)
		bad = numpy.where( resid >= o[-36])[0]
		
		from matplotlib import pyplot as plt
		fig = plt.figure()
		for i,b in enumerate(bad):
			l,m = dataDict['bls'][pol][b]
			print b, dataDict['bls'][pol][b], (aa.ants[l].stand, aa.ants[m].stand), resid[b]
			ax = fig.add_subplot(6, 6, i+1)
			ax.plot(numpy.angle(obsVis[b,:]) / (2*numpy.pi*fq), label='Data', marker='x')
			ax.plot(numpy.angle(simVis[b,:]) / (2*numpy.pi*fq), label='Model', marker='+')
			ax.plot(C[b,:] / (2*numpy.pi*fq), label='Angle', linestyle='--', marker='x')
			#ax.legend(loc=0)
		plt.draw()
		
		fig = plt.figure()
		ax = fig.gca()
		for i in xrange(C.shape[0]):
			ax.plot(C[i,:], marker='x')
		plt.draw()
		
		fig = plt.figure()
		for b in xrange(64):
			l,m = dataDict['bls'][pol][b]
			print b, dataDict['bls'][pol][b], (aa.ants[l].stand, aa.ants[m].stand), resid[b]
			ax = fig.add_subplot(8, 8, b+1)
			ax.plot(numpy.angle(obsVis[b,:]) / (2*numpy.pi*fq), label='Data', marker='x')
			ax.plot(numpy.angle(simVis[b,:]) / (2*numpy.pi*fq), label='Model', marker='+')
			ax.plot(C[b,:] / (2*numpy.pi*fq), label='Angle', linestyle='--', marker='x')
			#ax.legend(loc=0)
		plt.draw()
		
		plt.show()
	
	#C = C[:,len(chan)/2]
	fIndexList = range(0, fq.size-1, max([fq.size/N_CHAN,1]))
	if len(fIndexList) < 3:
		fIndexList = [0, fq.size/2, -1]
	
	Cp = numpy.zeros(len(fIndexList)*nBLs)
	for i,f in enumerate(fIndexList):
		for j in xrange(C.shape[0]):
			Cp[j+i*nBLs] = C[j,f]  / (2*numpy.pi*fq[f])
	
	return Cp


def selfCal(aa, dataDict, simDict, chan, pol, refAnt=0, returnDelays=False):
	"""Function used to perform a simple phase self-calibration of data stored in a 
	readUVData dictionary and a model sky stored in a lsl.sim.vis.buildSimSky 
	dictionary for a given polarization and channel(s)."""

	# Make sure we have the right polarization
	if pol not in dataDict['bls'].keys() and pol.lower() not in dataDict['bls'].keys():
		raise RuntimeError("Data dictionary does not have data for polarization '%s'" % pol)
	if pol not in simDict['bls'].keys() and pol.lower() not in simDict['bls'].keys():
		raise RuntimeError("Simulation dictionary does not have data for polarization '%s'" % pol)

	# Make sure that `chan' is an array by trying to find its length
	try:
		junk = len(chan)
	except TypeError:
		chan = [chan]

	N = len(aa.ants)
	found = False
	tauNames = []
	if refAnt != 0:
		origRefAnt = refAnt
		for i,ant in enumerate(aa.ants):
			if origRefAnt == ant.stand:
				print i, origRefAnt, refAnt, ant.stand
				refAnt = i
				found = True
			else:
				tauNames.append('tau%i' % ant.stand)
	else:
		found = True
	if not found:
		raise RuntimeError("Stand #%i not found in the array provided" % refAnt)
	print "Using antenna #%i as a reference (Stand #%i)" % (refAnt, aa.ants[refAnt].stand)
	
	nIter = 30
	
	tempGains = numpy.ones(N)
	tempDelays = numpy.zeros(N)
	
	histDelays = numpy.zeros((nIter,N))
	for i in xrange(nIter):
		##
		## Amplitude
		##
		#A = _buildAmplitudeA(aa, dataDict, simDict, chan, pol, refAnt=refAnt)
		#C = _buildAmplitudeC(aa, dataDict, simDict, chan, pol, refAnt=refAnt)
		#
		#good = numpy.where( numpy.isfinite(C) == 1 )[0]
		#A = A[good,:]
		#C = C[good]	
		#
		#bestGains, resid, rank, s = numpy.linalg.lstsq(A, C)
		#resid = numpy.array(C - numpy.dot(A, bestGains)).ravel()
		#resid = (C**2).sum(), (resid**2).sum()
		#bestGains = numpy.exp(bestGains)
		#tempGains *= bestGains
		
		## Report on the fits
		#print 'Best Gains: ', bestGains
		#print 'Residuals: ', resid
		#print 'Rank: ', rank
		
		#dataDict = scaleData(dataDict, bestGains, numpy.zeros_like(bestGains))
		
		#
		# Phase
		#
		
		A = _buildPhaseA(aa, dataDict, simDict, chan, pol, refAnt=refAnt)
		C = _buildPhaseC(aa, dataDict, simDict, chan, pol, refAnt=refAnt)
	
		good = numpy.where( numpy.isfinite(C) == 1 )[0]
		A = A[good,:]
		C = C[good]
		print C.shape, A.shape, len(tauNames), N
		
		m = ols.ols(numpy.array(C), numpy.array(A), y_varnm='deltaPhase', x_varnm=tauNames)
		m.summary()
		bestDelays = [m.b[j] for j in xrange(1, len(m.x_varnm))]
		
		#bestDelays, resid, rank, s = numpy.linalg.lstsq(A, C)
		#resid = numpy.array(C - numpy.dot(A, bestDelays)).ravel()
		#resid = (C**2).sum(), (resid**2).sum()
		
		bestDelays = list(bestDelays)
		bestDelays.insert(refAnt, m.b[0])
		bestDelays = numpy.array(bestDelays)
		tempDelays += bestDelays
		histDelays[i,:] = bestDelays
	
		# Report on the fits
		#print i, 'Best Delays: ', bestDelays
		print i, 'Range: ', bestDelays.min(), bestDelays.mean(), bestDelays.max()
		#print i, 'Residuals: ', resid
		#print i, 'Rank: ', rank
		
		dataDict = scaleData(dataDict, numpy.ones_like(bestDelays), bestDelays)
		if (numpy.abs(bestDelays)).max() <= 0.2:
			break
		
	#C = _buildPhaseC(aa, dataDict, simDict, chan, pol, refAnt=refAnt, plot=True)
	
	#import pylab
	#for i in xrange(50):
		#pylab.plot(histDelays[:,i].cumsum())
		#print i, m.x_varnm[i+1], tempDelays[i], (histDelays[:,i].cumsum())[-1]
	#pylab.show()
	
	print 'Best Gains: ', tempGains
	bestGains = tempGains
	
	print 'Best Delays: ', tempDelays
	bestDelays = tempDelays
	
	if returnDelays:
		return (dataDict, bestGains, bestDelays)
	else:
		return dataDict

