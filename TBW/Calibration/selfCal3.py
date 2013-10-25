# -*- coding: utf-8 -*-

"""Simple self-calibration module for correlated TBW and TBN data.

.. versionadded:: 0.5.5
"""

import sys
import math
import numpy

from lsl.imaging import utils
#from lsl.sim.vis import scaleData

__version__ = '0.1'
__revision__ = '$Rev: 497 $'
__all__ = ['selfCal', '__version__', '__revision__', '__all__']


def scaleData(dataDict, amps, delays, phases):
	"""
	Apply a set of antenna-based real gain values and phase delays in ns to a 
	data dictionary.  Returned the new scaled and delayed dictionary.
	
	..versionchanged:: 0.4.0
		The delays are now expected to be in nanoseconds rather than radians.
	"""

	import copy

	# Build the data dictionary to hold the scaled and delayed data
	sclUVData = {'freq': (dataDict['freq']).copy(), 'uvw': {}, 'vis': {}, 'wgt': {}, 'msk': {}, 'bls': {}, 'jd': {}}
	if dataDict['isMasked']:
		sclUVData['isMasked'] = True
	else:
		sclUVData['isMasked'] = False
	fq = dataDict['freq'] / 1e9
	
	cGains = []
	for i in xrange(len(amps)):
		cGains.append( amps[i]*numpy.exp(2j*numpy.pi*fq*delays[i] + 1j*phases[i]) )

	# Apply the scales and delays for all polarization pairs found in the original data
	for pol in dataDict['vis'].keys():
		sclUVData['bls'][pol] = []
		sclUVData['uvw'][pol] = []
		sclUVData['vis'][pol] = []
		sclUVData['wgt'][pol] = copy.copy(dataDict['wgt'][pol])
		sclUVData['msk'][pol] = copy.copy(dataDict['msk'][pol])
		sclUVData['jd'][pol] = copy.copy(dataDict['jd'][pol])

		for (i,j),uvw,vis in zip(dataDict['bls'][pol], dataDict['uvw'][pol], dataDict['vis'][pol]):
			sclUVData['bls'][pol].append( (i,j) )
			sclUVData['uvw'][pol].append( uvw )
			sclUVData['vis'][pol].append( vis*cGains[j].conj()*cGains[i] )

	return sclUVData


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
	fIndexList = range(0, fq.size-1, fq.size/100)
	if len(fIndexList) < 3:
		fIndexList = [0, fq.size/2, -1]
	
	A = numpy.zeros((len(fIndexList)*nBLs, 2*(nStands-1)))
	for i,f in enumerate(fIndexList):
		for j,(l,m) in enumerate(dataDict['bls'][pol]):
			if l < refAnt:
				A[j+i*nBLs,l]   =  2*numpy.pi*fq[f]
				A[j+i*nBLs,l + (nStands-1)] = 1 
			elif l > refAnt:
				A[j+i*nBLs,l-1] =  2*numpy.pi*fq[f]
				A[j+i*nBLs,l-1 + (nStands-1)] = 1
			else:
				pass
				
			if m < refAnt:
				A[j+i*nBLs,m]   = -2*numpy.pi*fq[f]
				A[j+i*nBLs,m + (nStands-1)]   = -1
			elif m > refAnt:
				A[j+i*nBLs,m-1] = -2*numpy.pi*fq[f]
				A[j+i*nBLs,m-1 + (nStands-1)] = -1
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
	
	C = numpy.angle(simVis / obsVis)
	
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
	fIndexList = range(0, fq.size-1, fq.size/100)
	if len(fIndexList) < 3:
		fIndexList = [0, fq.size/2, -1]
	
	Cp = numpy.zeros(len(fIndexList)*nBLs)
	for i,f in enumerate(fIndexList):
		for j in xrange(C.shape[0]):
			Cp[j+i*nBLs] = C[j,f]
	
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
	if refAnt != 0:
		for i,ant in enumerate(aa.ants):
			if refAnt == ant.stand:
				refAnt = i
				found = True
				break;
	else:
		found = True
	if not found:
		raise RuntimeError("Stand #%i not found in the array provided" % refAnt)
	print "Using antenna #%i as a reference (Stand #%i)" % (refAnt, aa.ants[refAnt].stand)

	tempGains = numpy.ones(N)
	tempDelays = numpy.zeros(N)
	tempPhases = numpy.zeros(N)
	for i in xrange(1):
		
		#
		# Amplitude
		#
		A = _buildAmplitudeA(aa, dataDict, simDict, chan, pol, refAnt=refAnt)
		C = _buildAmplitudeC(aa, dataDict, simDict, chan, pol, refAnt=refAnt)

		good = numpy.where( numpy.isfinite(C) == 1 )[0]
		A = A[good,:]
		C = C[good]	

		bestGains, resid, rank, s = numpy.linalg.lstsq(A, C)
		resid = numpy.array(C - numpy.dot(A, bestGains)).ravel()
		resid = (C**2).sum(), (resid**2).sum()
		bestGains = numpy.exp(bestGains)
		tempGains *= bestGains
		
		## Report on the fits
		#print 'Best Gains: ', bestGains
		#print 'Residuals: ', resid
		#print 'Rank: ', rank
		
		dataDict = scaleData(dataDict, bestGains, numpy.zeros_like(bestGains), numpy.zeros_like(bestGains))
		
		#
		# Phase
		#
		
		A = _buildPhaseA(aa, dataDict, simDict, chan, pol, refAnt=refAnt)
		C = _buildPhaseC(aa, dataDict, simDict, chan, pol, refAnt=refAnt)
	
		good = numpy.where( numpy.isfinite(C) == 1 )[0]
		A = A[good,:]
		C = C[good]
		
		bestDelays, resid, rank, s = numpy.linalg.lstsq(A, C)
		resid = numpy.array(C - numpy.dot(A, bestDelays)).ravel()
		resid = (C**2).sum(), (resid**2).sum()
		
		bestPhases = list(bestDelays)[len(bestDelays)/2:]
		bestDelays = list(bestDelays)[:len(bestDelays)/2]
		bestPhases.insert(refAnt, 0.0)
		bestDelays.insert(refAnt, 0.0)
		bestPhases = numpy.array(bestPhases)
		bestDelays = numpy.array(bestDelays)
		tempPhases += bestPhases
		tempDelays += bestDelays
	
		## Report on the fits
		#print 'Best Delays: ', bestDelays
		#print 'Residuals: ', resid
		#print 'Rank: ', rank
		
		dataDict = scaleData(dataDict, numpy.ones_like(bestDelays), bestDelays, bestPhases)
		
	#C = _buildPhaseC(aa, dataDict, simDict, chan, pol, refAnt=refAnt, plot=True)
	
	print 'Best Gains: ', tempGains
	bestGains = tempGains
	
	print 'Best Delays: ', tempDelays
	bestDelays = tempDelays
	
	print 'Best Phases: ', tempPhases
	bestPhases = tempPhases
	
	if returnDelays:
		return (dataDict, bestGains, bestDelays, bestPhases)
	else:
		return dataDict

