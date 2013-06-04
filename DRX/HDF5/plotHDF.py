#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Given a DRX HDF5 waterfall file, plot it in an interative way.

$Rev$
$LastChangedBy$
$LastChangedDate$
"""

import os
import sys
import h5py
import math
import time
import numpy
import ephem
import getopt
import subprocess
from datetime import datetime
from multiprocessing import Pool
from scipy.stats import scoreatpercentile as percentile

from lsl.common.dp import fS
from lsl.common import stations
from lsl.misc.mathutil import to_dB, from_dB, savitzky_golay
from lsl.statistics import robust
from lsl.statistics.kurtosis import spectralPower, std as skStd

import wx
import matplotlib
matplotlib.use('WXAgg')
matplotlib.interactive(True)

from matplotlib.backends.backend_wxagg import NavigationToolbar2WxAgg, FigureCanvasWxAgg
from matplotlib.figure import Figure


def usage(exitCode=None):
	print """plotHDF.py - Read in a HDF5 waterfall file and plot it interactively.

Usage: plotHDF.py [OPTIONS] file

Options:
-h, --help                  Display this help information
-s, --skip                  Skip the specified number of seconds at the beginning
                            of the file (default = 0)
-d, --duration              Number of seconds to display
                            (default = -1 to display all of the file)
-o, --observation           Plot the specified observation (default = 1)
"""

	if exitCode is not None:
		sys.exit(exitCode)
	else:
		return True


def parseOptions(args):
	config = {}
	# Command line flags - default values
	config['offset'] = 0.0
	config['duration'] = -1.0
	config['obsID'] = 1
	config['args'] = []

	# Read in and process the command line flags
	try:
		opts, args = getopt.getopt(args, "hs:d:o:", ["help", "skip=", "duration=", "observation="])
	except getopt.GetoptError, err:
		# Print help information and exit:
		print str(err) # will print something like "option -a not recognized"
		usage(exitCode=2)
	
	# Work through opts
	for opt, value in opts:
		if opt in ('-h', '--help'):
			usage(exitCode=0)
		elif opt in ('-s', '--skip'):
			config['offset'] = float(value)
		elif opt in ('-d', '--duration'):
			config['duration'] = float(value)
		elif opt in ('-o', '--observation'):
			config['obsID'] = int(value)
		else:
			assert False
	
	# Add in arguments
	config['args'] = args

	# Return configuration
	return config


def findMean(data):
	"""
	Tiny function to return the mean along the first axis.
	"""

	return numpy.mean(data, axis=0)


def findLimits(data, usedB=True):
	"""
	Tiny function to speed up the computing of the data range for the colorbar.
	Returns a two-element list of the lowest and highest values.
	"""
	
	dMin = data.min()
	if usedB:
		dMin = to_dB(dMin)
	if not numpy.isfinite(dMin):
		dMin = 0
		
	dMax = data.max()
	if usedB:
		dMax = to_dB(dMax)
	if not numpy.isfinite(dMax):
		dMax = dMin + 1
		
	return [dMin, dMax]


def bestFreqUnits(freq):
	"""Given a numpy array of frequencies in Hz, return a new array with the
	frequencies in the best units possible (kHz, MHz, etc.)."""

	# Figure out how large the data are
	try:
		scale = int(math.log10(freq.max()))
	except AttributeError:
		scale = int(math.log10(freq))
	if scale >= 9:
		divis = 1e9
		units = 'GHz'
	elif scale >= 6:
		divis = 1e6
		units = 'MHz'
	elif scale >= 3:
		divis = 1e3
		units = 'kHz'
	else:
		divis = 1
		units = 'Hz'

	# Convert the frequency
	newFreq = freq / divis

	# Return units and freq
	return (newFreq, units)


class Waterfall_GUI(object):
	def __init__(self, frame, freq=None, spec=None, tInt=None):
		self.frame = frame
		self.press = None
		
		self.filename = ''
		self.index = 0
		
		self.bandpass = False
		self.freq1 = freq
		self.freq2 = freq
		self.spec = spec
		self.tInt = None
		
		self.ax1a = None
		self.ax1b = None
		self.ax1c = None
		self.ax2 = None
		
		self.spectrumClick = None
		
		self.bandpassCut = 0.85
		self.driftOrder  = 5
		self.driftCut    = 4
		self.kurtosisSec = 1
		self.kurtosisCut = 4
		
	def loadData(self, filename, obsID=1):
		"""
		Load in data from an NPZ file.
		"""
		
		print "Loading file '%s'" % os.path.split(filename)[1]
		tStart = time.time()
		
		# Save the filename
		self.filename = filename
		self.obsID = obsID
		h = h5py.File(self.filename, 'r')
		obs = h.get('Observation%i' % obsID, None)
		if obs is None:
			raise RuntimeError("No observation #%i in file '%s'" % (obsID, os.path.basename(self.filename)))
		
		# Load the Data
		print " %6.3f s - Extracting data for observation #%i" % (time.time() - tStart, obsID)
		self.beam  = obs.attrs['Beam']
		self.srate = obs.attrs['sampleRate']
		self.tInt  = obs.attrs['tInt']
		self.time  = numpy.zeros(obs['time'].shape, dtype=obs['time'].dtype)
		obs['time'].read_direct(self.time)
		self.time -= self.time[0]
		valid = numpy.where( self.time >= 0 )
		self.time = self.time[valid]
		
		tuning1 = obs.get('Tuning1', None)
		tuning2 = obs.get('Tuning2', None)
		
		dataProducts = list(tuning1)
		mapper = {'XX': 0, 'I': 0, 'XY': 1, 'Q': 1, 'YX': 2, 'U': 2, 'YY': 3, 'V': 3}
		for exclude in ('freq', 'Mask', 'Saturation', 'SpectralKurtosis'):
			try:
				ind = dataProducts.index(exclude)
				del dataProducts[ind]
			except ValueError:
				pass
		if dataProducts[0][0] in ('X', 'Y'):
			self.linear = True
			self.usedB = True
		else:
			self.linear = False
			self.usedB = False
			
		# Figure out the selection process
		self.iOffset = int(round(self.frame.offset / self.tInt))
		if self.frame.duration < 0:
			self.iDuration = self.time.size - self.iOffset
		else:
			self.iDuration = int(round(self.frame.duration / self.tInt))
		## Make sure we don't fall off the end of the file
		if self.iOffset + self.iDuration > tuning1[dataProducts[0]].shape[0]:
			self.iDuration = tuning1[dataProducts[0]].shape[0] - self.iOffset
		selection = numpy.s_[self.iOffset:self.iOffset+self.iDuration, :]
		
		if self.iOffset != 0:
			print "            -> Offsetting %i integrations (%.3f s)" % (self.iOffset, self.iOffset*self.tInt)
		print "            -> Displaying %i integrations (%.3f s)" % (self.iDuration, self.iDuration*self.tInt)
		
		self.time = self.time[:self.iDuration]
		
		self.freq1 = numpy.zeros(tuning1['freq'].shape, dtype=tuning1['freq'].dtype)
		tuning1['freq'].read_direct(self.freq1)
		self.freq2 = numpy.zeros(tuning2['freq'].shape, dtype=tuning2['freq'].dtype)
		tuning2['freq'].read_direct(self.freq2)
		
		self.spec = numpy.empty((self.iDuration, 8, self.freq1.size), dtype=numpy.float32)
		
		for p in dataProducts:
			## Tuning 1
			ind = 4*(1-1) + mapper[p]
			part = numpy.empty((self.iDuration, self.freq1.size), dtype=tuning1[p].dtype)
			tuning1[p].read_direct(part, selection)
			self.spec[:,ind,:] = part.astype(numpy.float32)
			
			## Tuning 2
			ind = 4*(2-1) + mapper[p]
			part = numpy.empty((self.iDuration, self.freq2.size), dtype=tuning2[p].dtype)
			tuning2[p].read_direct(part, selection)
			self.spec[:,ind,:] = part.astype(numpy.float32)
			
			del part
			
		self.sats = numpy.empty((self.iDuration, 4), dtype=numpy.float32)
		part = numpy.empty((self.iDuration, 2), dtype=tuning1['Saturation'].dtype)
		tuning1['Saturation'].read_direct(part, selection)
		self.sats[:,0:2] = part / (self.tInt * self.srate)
		part = numpy.empty((self.iDuration, 2), dtype=tuning2['Saturation'].dtype)
		tuning2['Saturation'].read_direct(part, selection)
		self.sats[:,2:4] = part / (self.tInt * self.srate)
		del part
		# Mask out negative saturation values since that indicates the data is
		# not avaliable
		self.sats = numpy.ma.array(self.sats, mask=(self.sats < 0))
		
		mask1 = tuning1.get('Mask', None)
		mask2 = tuning2.get('Mask', None)
		
		mask = numpy.zeros(self.spec.shape, dtype=numpy.bool)
		
		for p in dataProducts:
			if mask1 is not None:
				## Tuning 1
				ind = 4*(1-1) + mapper[p]
				part = numpy.empty((self.iDuration, self.freq1.size), dtype=mask1[p].dtype)
				mask1[p].read_direct(part, selection)
				mask[:,ind,:] = part.astype(numpy.float32)
				
				del part
				
			if mask2 is not None:
				## Tuning 2
				ind = 4*(2-1) + mapper[p]
				part = numpy.empty((self.iDuration, self.freq2.size), dtype=mask2[p].dtype)
				mask2[p].read_direct(part, selection)
				mask[:,ind,:] = part.astype(numpy.float32)
				
				del part
		
		self.spec = numpy.ma.array(self.spec, mask=mask)
		
		# Construct frequency and time master masks to prevent some masked things from getting unmasked
		self.freqMask = numpy.median(self.spec.mask, axis=0)
		self.timeMask = numpy.median(self.spec.mask, axis=2)
		
		# Other data to keep around
		self.timesNPZ = numpy.zeros(obs['time'].shape, dtype=obs['time'].dtype)
		obs['time'].read_direct(self.timesNPZ)
		
		# Deal with the potential for aggregated files
		self.tIntActual = self.tInt
		self.tIntOriginal = self.tInt
		self.filenames = None
		
		# Gather up the target information
		self.dataProducts = dataProducts
		self.target = obs.attrs['TargetName']
		self.ra = obs.attrs['RA']
		self.raUnit = obs.attrs['RA_Units']
		self.dec = obs.attrs['Dec']
		self.decUnit = obs.attrs['Dec_Units']
		self.mode = obs.attrs['TrackingMode']
		self.rbw = obs.attrs['RBW']
		self.rbwUnit = obs.attrs['RBW_Units']
		
		# Close out the file
		h.close()
		
		## Get the filter model
		#print " %6.3f s - Building DRX bandpass model" % (time.time() - tStart)
		#self.bpm = drxFilter(sampleRate=self.srate)(self.freq1)
		
		# Compute the bandpass fit
		print " %6.3f s - Computing bandpass fits" % (time.time() - tStart)
		self.computeBandpass()
		
		# Find the mean spectra
		print " %6.3f s - Computing mean spectra" % (time.time() - tStart)
		try:
			from _helper import FastAxis0Mean
			self.mean = FastAxis0Mean(self.spec)
			self.meanBandpass = FastAxis0Mean(self.specBandpass)
		except ImportError:
			self.mean = numpy.mean(self.spec, axis=0)
			self.meanBandpass = numpy.mean(self.specBandpass, axis=0)
		
		# Set default colobars
		print " %6.3f s - Setting default colorbar ranges" % (time.time() - tStart)
		self.limits = [None,]*self.spec.shape[1]
		self.limitsBandpass = [None,]*self.spec.shape[1]
		
		try:
			from _helper import FastAxis1MinMax
			limits0 = FastAxis1MinMax(self.spec)
			limits1 = FastAxis1MinMax(self.specBandpass, chanMin=self.spec.shape[2]/10, chanMax=9*self.spec.shape[2]/10)
			if self.usedB:
				limits0 = to_dB(limits0)
				limits1 = to_dB(limits1)
			for i in xrange(self.spec.shape[1]):
				self.limits[i] = list(limits0[i,:])
				self.limitsBandpass[i] = list(limits1[i,:])
		except ImportError:
			for i in xrange(self.spec.shape[1]):
				self.limits[i] = findLimits(self.spec[:,i,:], usedB=self.usedB)
			for i in xrange(self.spec.shape[1]):
				self.limitsBandpass[i] = findLimits(self.specBandpass[:,i,toUse], usedB=self.usedB)
				
		try:
			self.disconnect()
		except:
			pass
			
		print " %6.3f s - Finished preparing data" % (time.time() - tStart)
		
	def render(self):
		# Clear the old marks
		self.oldMarkA = None
		self.oldMarkB = None
		self.oldMarkC = None
		
		# Clear the old figures
		self.frame.figure1a.clf()
		self.frame.figure1b.clf()
		self.frame.figure1c.clf()
		self.frame.figure2.clf()
		
		self.connect()
		
	def computeBandpass(self):
		"""
		Compute the bandpass fits.
		"""
		
		try:
			from _helper import FastAxis0Median
			meanSpec = FastAxis0Median(self.spec.data)
		except ImportError:
			meanSpec = numpy.median(self.spec.data, axis=0)
			
		bpm2 = []
		for i in xrange(self.spec.shape[1]):
			bpm = savitzky_golay(meanSpec[i,:], 41, 9, deriv=0)
			
			if bpm.mean() == 0:
				bpm += 1
			bpm2.append( bpm / bpm.mean() )
			
		# Apply the bandpass correction
		bpm2 = numpy.array(bpm2)
		self.specBandpass = numpy.ma.array(self.spec.data*1.0, mask=self.spec.mask)
		try:
			from _helper import FastAxis0Bandpass
			FastAxis0Bandpass(self.specBandpass.data, bpm2.astype(numpy.float32))
		except ImportError:
			for i in xrange(self.spec.shape[1]):
				self.specBandpass[:,i,:] = self.spec[:,i,j] / bpm2[i][j]
				
		return True
		
	def draw(self):
		"""
		Draw the waterfall diagram and the total power with time.
		"""
		
		if self.index / (self.spec.shape[1]/2) == 0:
			freq = self.freq1
		else:
			freq = self.freq2
		
		if self.bandpass:
			spec = self.specBandpass
			limits = self.limitsBandpass
		else:
			spec = self.spec
			limits = self.limits
		
		# Plot 1(a) - Waterfall
		self.frame.figure1a.clf()
		self.ax1a = self.frame.figure1a.gca()
		if self.usedB:
			m = self.ax1a.imshow(to_dB(spec[:,self.index,:]), interpolation='nearest', extent=(freq[0]/1e6, freq[-1]/1e6, self.time[0], self.time[-1]), origin='lower', vmin=limits[self.index][0], vmax=limits[self.index][1])
			try:
				cm = self.frame.figure1a.colorbar(m, use_gridspec=True)
			except:
				cm = self.frame.figure1a.colorbar(m)
			cm.ax.set_ylabel('PSD [arb. dB]')
		else:
			m = self.ax1a.imshow(spec[:,self.index,:], interpolation='nearest', extent=(freq[0]/1e6, freq[-1]/1e6, self.time[0], self.time[-1]), origin='lower', vmin=limits[self.index][0], vmax=limits[self.index][1])
			try:
				cm = self.frame.figure1a.colorbar(m, use_gridspec=True)
			except:
				cm = self.frame.figure1a.colorbar(m)
			cm.ax.set_ylabel('PSD [arb. lin.]')
		self.ax1a.axis('auto')
		self.ax1a.set_xlim((freq[0]/1e6, freq[-1]/1e6))
		self.ax1a.set_ylim((self.time[0], self.time[-1]))
		self.ax1a.set_xlabel('Frequency [MHz]')
		self.ax1a.set_ylabel('Elapsed Time - %.3f [s]' % (self.iOffset*self.tInt))
		if self.linear:
			tun = self.index / 4 + 1
			ind = self.index % 4
			mapper = {0: 'XX', 1: 'XY', 2: 'YX', 3: 'YY'}
			self.ax1a.set_title('Tuning %i, %s' % (tun, mapper[ind]))
		else:
			tun = self.index / 4 + 1
			ind = self.index % 4
			mapper = {0: 'I', 1: 'Q', 2: 'U', 3: 'V'}
			self.ax1a.set_title('Tuning %i, %s' % (tun, mapper[ind]))
			
		if self.oldMarkA is not None:
			self.ax1a.lines.extend(self.oldMarkA)
		
		try:
			self.frame.figure1a.tight_layout()
		except:
			pass
		self.frame.canvas1a.draw()
		
		# Plot 1(b) - Saturation Fraction
		self.frame.figure1b.clf()
		self.ax1b = self.frame.figure1b.gca()
		self.ax1b.plot(self.sats[:,2*(tun-1)+0], self.time, linestyle='-', color='blue')
		self.ax1b.plot(self.sats[:,2*(tun-1)+1], self.time, linestyle='-', color='green')
		self.ax1b.set_xlim((-0.05, 1.05))
		self.ax1b.set_ylim((self.time[0], self.time[-1]))
		self.ax1b.set_xlabel('Saturation Fraction')
		self.ax1b.set_ylabel('Elapsed Time - %.3f [s]' % (self.iOffset*self.tInt))
		self.ax1b.xaxis.set_ticks([0.0, 0.25, 0.5, 0.75, 1.0])
		self.ax1b.xaxis.set_ticklabels(['0', '', '0.5', '', '1'])
		
		if self.oldMarkB is not None:
			self.ax1b.lines.extend(self.oldMarkB)
		try:
			self.frame.figure1b.tight_layout()
		except:
			pass
		self.frame.canvas1b.draw()
		
		# Plot 1(c) - Drift
		self.drift = spec[:,:,spec.shape[2]/8:7*spec.shape[2]/8].sum(axis=2)
		
		self.frame.figure1c.clf()
		self.ax1c = self.frame.figure1c.gca()
		if self.usedB:
			self.ax1c.plot(to_dB(self.drift[:,self.index]), self.time, linestyle=' ', marker='x')
			self.ax1c.set_xlabel('Inner 75% Total Power [arb. dB]')
		else:
			self.ax1c.plot(self.drift[:,self.index], self.time, linestyle=' ', marker='x')
			self.ax1c.set_xlabel('Inner 75% Total Power [arb. lin.]')
		self.ax1c.set_ylim((self.time[0], self.time[-1]))
		self.ax1c.set_ylabel('Elapsed Time - %.3f [s]' % (self.iOffset*self.tInt))
		
		if self.oldMarkC is not None:
			self.ax1c.lines.extend(self.oldMarkC)
		
		try:
			self.frame.figure1c.tight_layout()
		except:
			pass
		self.frame.canvas1c.draw()
		
	def drawSpectrum(self, clickY):
		"""Get the spectrum at a particular point in time."""
		
		try:
			dataY = int(round(clickY / self.tInt))
		except TypeError:
			return False
			
		if self.index / (self.spec.shape[1]/2) == 0:
			freq = self.freq1
		else:
			freq = self.freq2
		
		if self.bandpass:
			spec = self.specBandpass[dataY,self.index,:]
			medianSpec = self.meanBandpass[self.index,:]
			limits = self.limitsBandpass
		else:
			spec = self.spec[dataY,self.index,:]
			medianSpec = self.mean[self.index,:]
			limits = self.limits
		
		if self.frame.toolbar.mode == 'zoom rect':
			try:
				oldXlim = self.ax2.get_xlim()
				oldYlim = self.ax2.get_ylim()
			except:
				oldXlim = [freq[0]/1e6, freq[-1]/1e6]
				oldYlim = limits[self.index]
		else:
			oldXlim = [freq[0]/1e6, freq[-1]/1e6]
			oldYlim = limits[self.index]
		
		self.frame.figure2.clf()
		self.ax2 = self.frame.figure2.gca()
		if self.usedB:
			self.ax2.plot(freq/1e6, to_dB(spec), linestyle=' ', marker='o', label='Current', mec='blue', mfc='None')
			self.ax2.plot(freq/1e6, to_dB(medianSpec), label='Mean', alpha=0.5, color='green')
			self.ax2.set_ylabel('PSD [arb. dB]')
		else:
			self.ax2.plot(freq/1e6, spec, linestyle=' ', marker='o', label='Current', mec='blue', mfc='None')
			self.ax2.plot(freq/1e6, medianSpec, label='Mean', alpha=0.5, color='green')
			self.ax2.set_ylabel('PSD [arb. lin.]')
		self.ax2.set_xlim(oldXlim)
		self.ax2.set_ylim(oldYlim)
		self.ax2.legend(loc=0)
		self.ax2.set_xlabel('Frequency [MHz]')
		
		if self.filenames is None:
			if self.bandpass:
				self.ax2.set_title("%s UTC + bandpass" % datetime.utcfromtimestamp(self.timesNPZ[dataY]))
			else:
				self.ax2.set_title("%s UTC" % datetime.utcfromtimestamp(self.timesNPZ[dataY]))
		else:
			if self.bandpass:
				self.ax2.set_title("%s + bandpass" % self.filenames[dataY])
			else:
				self.ax2.set_title(self.filenames[dataY])
		
		try:
			self.frame.figure2.tight_layout()
		except:
			pass
		self.frame.canvas2.draw()
		self.spectrumClick = clickY
		
	def makeMark(self, clickY):
		
		try:
			dataY = int(round(clickY / self.tInt))
		except TypeError:
			return False
		
		if self.oldMarkA is not None:
			try:
				del self.ax1a.lines[-1]
			except:
				pass
			
		if self.oldMarkB is not None:
			try:
				del self.ax1b.lines[-1]
			except:
				pass
				
		if self.oldMarkC is not None:
			try:
				del self.ax1c.lines[-1]
			except:
				pass
		
		oldXSizeA = self.ax1a.get_xlim()
		oldYSizeA = self.ax1a.get_ylim()
		
		oldXSizeB = self.ax1b.get_xlim()
		oldYSizeB = self.ax1b.get_ylim()
		
		oldXSizeC = self.ax1c.get_xlim()
		oldYSizeC = self.ax1c.get_ylim()
		
		self.oldMarkA = self.ax1a.plot([-1e20, 1e20], [self.time[dataY],]*2, color='red')
		self.oldMarkB = self.ax1b.plot([-1e20, 1e20], [self.time[dataY],]*2, color='red')
		self.oldMarkC = self.ax1c.plot([-1e20, 1e20], [self.time[dataY],]*2, color='red')
		
		self.ax1a.set_xlim(oldXSizeA)
		self.ax1a.set_ylim(oldYSizeA)
		
		self.ax1b.set_xlim(oldXSizeB)
		self.ax1b.set_ylim(oldYSizeB)
		
		self.ax1c.set_xlim(oldXSizeC)
		self.ax1c.set_ylim(oldYSizeC)
		
		self.frame.canvas1a.draw()
		self.frame.canvas1b.draw()
		self.frame.canvas1c.draw()
		
	def suggestMask(self, index):
		"""
		Suggest a mask for the current index.
		"""
		
		if self.index / (self.spec.shape[1]/2) == 0:
			freq = self.freq1
		else:
			freq = self.freq2
		freqP = freq - freq[len(freq)/2]
		
		bad = numpy.where( (freqP < -self.srate/2*self.bandpassCut) | (freqP > self.srate/2*self.bandpassCut) )[0]
		for b in bad:
			self.spec.mask[:,index,b] = True
			self.specBandpass.mask[:,index,b] = True
			self.freqMask[index,b] = True
		
		spec = self.spec.data[:,index,:]
		drift = spec.sum(axis=1)
		coeff = robust.polyfit(numpy.arange(drift.size), drift, self.driftOrder)
		fit = numpy.polyval(coeff, numpy.arange(drift.size))
		rDrift = drift / fit
		
		mean = robust.mean(rDrift)
		std  = robust.std(rDrift)
		
		bad = numpy.where( numpy.abs(rDrift - mean) >= self.driftCut*std )[0]
		for b in bad:
			self.spec.mask[b,index,:] = True
			self.specBandpass.mask[b,index,:] = True
			self.timeMask[b,index] = True
		
		N = self.srate/(freq.size+1)*self.tIntActual
		kurtosis = numpy.zeros((self.kurtosisSec, self.spec.shape[2]))
		
		secSize = self.spec.shape[0]/self.kurtosisSec
		for k in xrange(self.kurtosisSec):
			tStart = k*secSize
			tStop  = (k+1)*secSize
			
			for j in xrange(self.spec.shape[2]):
				channel = self.spec.data[tStart:tStop,index,j]
				kurtosis[k,j] = spectralPower(channel, N=N)
		
		kMean = 1.0
		kStd  = skStd(secSize, N)
		# Correction for some averaged data sets - probably not all that vaid
		# for more robust things.  However, this only hits the aggregated files
		# since `self.filenames` is None for stand-alone files.
		if self.filenames is not None:
			kurtosis /= kurtosis.mean()
		
		bad = numpy.where( numpy.abs(kurtosis - kMean) >= self.kurtosisCut*kStd )
		for k,b in zip(bad[0], bad[1]):
			tStart = k*secSize
			tStop  = (k+1)*secSize
			
			try:
				for j in xrange(b-2, b+3):
					self.spec.mask[tStart:tStop,index,j] = True
					self.specBandpass.mask[tStart:tStop,index,j] = True
					self.freqMask[index,j] = True
			except IndexError:
				pass
			
		return True
		
	def resetMask(self, index):
		"""
		Reset the specified mask.
		"""
		
		self.spec.mask[:,index,:] = False
		self.specBandpass.mask[:,index,:]= False
		self.timeMask[:,index] = False
		self.freqMask[index,:] = False
		
		return True
		
	def connect(self):
		"""
		Connect to all the events we need
		"""
		
		self.cidpress1a = self.frame.figure1a.canvas.mpl_connect('button_press_event', self.on_press1a)
		self.cidpress1b = self.frame.figure1b.canvas.mpl_connect('button_press_event', self.on_press1b)
		self.cidpress1c = self.frame.figure1c.canvas.mpl_connect('button_press_event', self.on_press1c)
		self.cidpress2  = self.frame.figure2.canvas.mpl_connect('button_press_event', self.on_press2)
		self.cidmotion  = self.frame.figure1a.canvas.mpl_connect('motion_notify_event', self.on_motion)
		self.frame.toolbar.update = self.on_update
		
	def on_update(self, *args):
		"""
		Override the toolbar 'update' operation.
		"""
		
		self.frame.toolbar.set_history_buttons()
		
	def on_press1a(self, event):
		"""
		On button press we will see if the mouse is over us and store some data
		"""
		
		if event.inaxes:
			clickX = event.xdata
			clickY = event.ydata
			
			dataY = int(round(clickY / self.tInt))
			
			if event.button == 1:
				self.drawSpectrum(clickY)
				self.makeMark(clickY)
			elif event.button == 2:
				self.spec.mask[dataY, self.index, :] = self.freqMask[self.index,:]
				self.specBandpass.mask[dataY, self.index, :] = self.freqMask[self.index,:]
				self.timeMask[dataY, self.index] = False
				self.draw()
			elif event.button == 3:
				self.spec.mask[dataY, self.index, :] = True
				self.specBandpass.mask[dataY, self.index, :] = True
				self.timeMask[dataY, self.index] = True
				self.draw()
			else:
				pass
				
	def on_press1b(self, event):
		"""
		On button press we will see if the mouse is over us and store some data
		"""
		
		if event.inaxes:
			clickX = event.xdata
			clickY = event.ydata
			
			dataY = int(round(clickY / self.tInt))
			
			if event.button == 1:
				self.drawSpectrum(clickY)
				self.makeMark(clickY)
			elif event.button == 2:
				self.spec.mask[dataY, self.index, :] = self.freqMask[self.index,:]
				self.specBandpass.mask[dataY, self.index, :] = self.freqMask[self.index,:]
				self.timeMask[dataY, self.index] = False
				self.draw()
			elif event.button == 3:
				self.spec.mask[dataY, self.index, :] = True
				self.specBandpass.mask[dataY, self.index, :] = True
				self.timeMask[dataY, self.index] = True
				self.draw()
			else:
				pass
				
	def on_press1c(self, event):
		"""
		On button press we will see if the mouse is over us and store some data
		"""
		
		if event.inaxes:
			clickX = event.xdata
			clickY = event.ydata
			
			scaleX = self.ax1c.get_xlim()
			rangeX = scaleX[1] - scaleX[0]
			
			scaleY = self.ax1c.get_ylim()
			rangeY = scaleY[1] - scaleY[0]
			
			dataY = int(round(clickY / self.tInt))
			lower = dataY - 200
			lower = 0 if lower < 0 else lower
			upper = dataY + 200
			upper = self.drift.shape[0]-1 if upper > self.drift.shape[0]-1 else upper
			
			if self.usedB:
				d =  ((clickX - to_dB(self.drift.data[lower:upper,self.index]))/rangeX)**2
				d += ((clickY - self.time[lower:upper])/rangeY)**2
				d = numpy.sqrt(d)
				best = numpy.where( d == d.min() )[0][0] + lower
				bestD = d[best - lower]
				
				print "Clicked at %.3f, %.3f => resolved to entry %i at %.3f, %.3f" % (clickX, clickY, best, to_dB(self.drift.data[best, self.index]), self.time[best])
			else:
				d =  ((clickX - self.drift.data[lower:upper,self.index])/rangeX)**2
				d += ((clickY - self.time[lower:upper])/rangeY)**2
				d = numpy.sqrt(d)
				best = numpy.where( d == d.min() )[0][0] + lower
				bestD = d[best - lower]
				
				print "Clicked at %.3f, %.3f => resolved to entry %i at %.3f, %.3f" % (clickX, clickY, best, self.drift.data[best, self.index], self.time[best])
				
			if event.button == 1:
				self.drawSpectrum(clickY)
				self.makeMark(clickY)
			elif event.button == 2:
				self.spec.mask[best, self.index, :] = self.freqMask[self.index,:]
				self.specBandpass.mask[best, self.index, :] = self.freqMask[self.index,:]
				self.timeMask[best, self.index] = False
				self.draw()
			elif event.button == 3:
				self.spec.mask[best, self.index, :] = True
				self.specBandpass.mask[best, self.index, :] = True
				self.timeMask[best, self.index] = True
				self.draw()
			else:
				pass
			
	def on_press2(self, event):
		"""
		On button press we will see if the mouse is over us and store some data
		"""
		
		if event.inaxes:
			clickX = event.xdata
			clickY = event.ydata
			
			if self.index / (self.spec.shape[1]/2) == 0:
				freq = self.freq1
			else:
				freq = self.freq2
			
			dataX = numpy.where(numpy.abs(clickX-freq/1e6) == (numpy.abs(clickX-freq/1e6).min()))[0][0]
			
			if event.button == 2:
				self.spec.mask[:, self.index, dataX] = self.timeMask[:,self.index]
				self.specBandpass.mask[:, self.index, dataX] = self.timeMask[:,self.index]
				self.freqMask[self.index, dataX] = False
			elif event.button == 3:
				self.spec.mask[:, self.index, dataX] = True
				self.specBandpass.mask[:, self.index, dataX] = True
				self.freqMask[self.index, dataX] = True
			else:
				pass
			
			self.draw()
			self.drawSpectrum(self.spectrumClick)
			
	def on_motion(self, event):
		"""
		On mouse motion display the data value under the cursor
		"""
		
		if event.inaxes:
			clickX = event.xdata
			clickY = event.ydata
			
			if self.index / (self.spec.shape[1]/2) == 0:
				freq = self.freq1
			else:
				freq = self.freq2
			
			dataX = numpy.where(numpy.abs(clickX-freq/1e6) == (numpy.abs(clickX-freq/1e6).min()))[0][0]
			dataY = numpy.where(numpy.abs(clickY-self.time) == (numpy.abs(clickY-self.time).min()))[0][0]
			if not self.spec.mask[dataY, self.index, dataX]:
				if self.usedB:
					value = to_dB(self.spec[dataY, self.index, dataX])
					self.frame.statusbar.SetStatusText("f=%.4f MHz, t=%.4f s, p=%.2f dB" % (clickX, clickY, value))
				else:
					value = self.spec[dataY, self.index, dataX]
					self.frame.statusbar.SetStatusText("f=%.4f MHz, t=%.4f s, p=%.2f" % (clickX, clickY, value))
			else:
				self.frame.statusbar.SetStatusText("")
		else:
			self.frame.statusbar.SetStatusText("")
			
	def disconnect(self):
		"""
		Disconnect all the stored connection ids
		"""
		
		self.frame.figure1a.canvas.mpl_disconnect(self.cidpress1a)
		self.frame.figure1b.canvas.mpl_disconnect(self.cidpress1b)
		self.frame.figure1c.canvas.mpl_disconnect(self.cidpress1c)
		self.frame.figure2.canvas.mpl_disconnect(self.cidpress2)
		self.frame.figure1a.canvas.mpl_disconnect(self.cidmotion)


ID_OPEN    = 10
ID_SAVE    = 11
ID_SAVE_AS = 12
ID_QUIT    = 13

ID_COLOR_AUTO = 20
ID_COLOR_ADJUST = 21

ID_TUNING1_1 = 30
ID_TUNING1_2 = 31
ID_TUNING1_3 = 32
ID_TUNING1_4 = 33
ID_TUNING2_1 = 34
ID_TUNING2_2 = 35
ID_TUNING2_3 = 36
ID_TUNING2_4 = 37
ID_CHANGE_RANGE = 38
ID_CHANGE_OBSID = 39

ID_MASK_SUGGEST_CURRENT = 40
ID_MASK_SUGGEST_ALL = 41
ID_MASK_RESET_CURRENT = 42
ID_MASK_RESET_ALL = 43
ID_MASK_TWEAK = 44

ID_BANDPASS_ON = 50
ID_BANDPASS_OFF = 51
ID_BANDPASS_RECOMPUTE = 52

ID_DETAIL_SUMMARY = 60
ID_DETAIL_SUBFILE = 61
ID_DETAIL_WATERFALL = 62
ID_DETAIL_DRIFTCURVE = 63

class MainWindow(wx.Frame):
	def __init__(self, parent, id):
		self.dirname = ''
		self.filename = ''
		self.offset = 0.0
		self.duration = -1
		self.data = None
		
		self.examineWindow = None
		
		wx.Frame.__init__(self, parent, id, title="DRX Waterfall Viewer", size=(1000,600))
		
	def render(self):
		self.initUI()
		self.initEvents()
		self.Show()
		self.SetClientSize((1000,600))
		
		self.cAdjust = None
		
	def initUI(self):
		self.statusbar = self.CreateStatusBar() # A Statusbar in the bottom of the window
		
		font = wx.SystemSettings_GetFont(wx.SYS_SYSTEM_FONT)
		font.SetPointSize(10)
		
		menuBar = wx.MenuBar()
		
		fileMenu = wx.Menu()
		colorMenu = wx.Menu()
		dataMenu = wx.Menu()
		maskMenu = wx.Menu()
		bandpassMenu = wx.Menu()
		detailsMenu = wx.Menu()
		
		## File Menu
		open = wx.MenuItem(fileMenu, ID_OPEN, "&Open")
		#fileMenu.AppendItem(open)
		save = wx.MenuItem(fileMenu, ID_SAVE, "&Save")
		fileMenu.AppendItem(save)
		saveas = wx.MenuItem(fileMenu, ID_SAVE_AS, "Save &As")
		fileMenu.AppendItem(saveas)
		fileMenu.AppendSeparator()
		exit = wx.MenuItem(fileMenu, ID_QUIT, "E&xit")
		fileMenu.AppendItem(exit)
		
		## Color Menu
		auto = wx.MenuItem(colorMenu, ID_COLOR_AUTO, '&Auto-scale Colorbar')
		colorMenu.AppendItem(auto)
		cadj = wx.MenuItem(colorMenu, ID_COLOR_ADJUST, '&Adjust Contrast')
		colorMenu.AppendItem(cadj)
		
		## Data Menu
		if self.data.linear:
			t1p1 = dataMenu.AppendRadioItem(ID_TUNING1_1, 'Tuning 1, XX')
			t1p2 = dataMenu.AppendRadioItem(ID_TUNING1_2, 'Tuning 1, XY')
			t1p3 = dataMenu.AppendRadioItem(ID_TUNING1_3, 'Tuning 1, YX')
			t1p4 = dataMenu.AppendRadioItem(ID_TUNING1_4, 'Tuning 1, YY')
			t2p1 = dataMenu.AppendRadioItem(ID_TUNING2_1, 'Tuning 2, XX')
			t2p2 = dataMenu.AppendRadioItem(ID_TUNING2_2, 'Tuning 2, XY')
			t2p3 = dataMenu.AppendRadioItem(ID_TUNING2_3, 'Tuning 2, YX')
			t2p4 = dataMenu.AppendRadioItem(ID_TUNING2_4, 'Tuning 2, YY')
		else:
			t1p1 = dataMenu.AppendRadioItem(ID_TUNING1_1, 'Tuning 1, I')
			t1p2 = dataMenu.AppendRadioItem(ID_TUNING1_2, 'Tuning 1, Q')
			t1p3 = dataMenu.AppendRadioItem(ID_TUNING1_3, 'Tuning 1, U')
			t1p4 = dataMenu.AppendRadioItem(ID_TUNING1_4, 'Tuning 1, V')
			t2p1 = dataMenu.AppendRadioItem(ID_TUNING2_1, 'Tuning 2, I')
			t2p2 = dataMenu.AppendRadioItem(ID_TUNING2_2, 'Tuning 2, Q')
			t2p3 = dataMenu.AppendRadioItem(ID_TUNING2_3, 'Tuning 2, U')
			t2p4 = dataMenu.AppendRadioItem(ID_TUNING2_4, 'Tuning 2, V')
			
		dataMenu.InsertSeparator(4)
		dataMenu.AppendSeparator()
		self.changeRangeButton = wx.MenuItem(colorMenu, ID_CHANGE_RANGE, 'Change Time &Range')
		self.changeObservationButton = wx.MenuItem(colorMenu, ID_CHANGE_OBSID, 'Change &Observation')
		dataMenu.AppendItem(self.changeRangeButton)
		dataMenu.AppendItem(self.changeObservationButton)
		
		t1p1.Enable(False)
		t1p2.Enable(False)
		t1p3.Enable(False)
		t1p4.Enable(False)
		t2p1.Enable(False)
		t2p2.Enable(False)
		t2p3.Enable(False)
		t2p4.Enable(False)
		for p in self.data.dataProducts:
			if p in ('I', 'XX'):
				t1p1.Enable(True)
				t2p1.Enable(True)
			elif p in ('Q', 'XY'):
				t1p2.Enable(True)
				t2p2.Enable(True)
			elif p in ('U', 'YX'):
				t1p3.Enable(True)
				t2p3.Enable(True)
			else:
				t1p4.Enable(True)
				t2p4.Enable(True)
				
		## Mask Menu
		suggestC = wx.MenuItem(maskMenu, ID_MASK_SUGGEST_CURRENT, 'Suggest Mask - Current')
		maskMenu.AppendItem(suggestC)
		suggestA = wx.MenuItem(maskMenu, ID_MASK_SUGGEST_ALL, 'Suggest Mask - All')
		maskMenu.AppendItem(suggestA)
		maskMenu.AppendSeparator()
		resetC = wx.MenuItem(maskMenu, ID_MASK_RESET_CURRENT, 'Reset Mask - Current')
		maskMenu.AppendItem(resetC)
		resetA = wx.MenuItem(maskMenu, ID_MASK_RESET_ALL, 'Reset Mask - All')
		maskMenu.AppendItem(resetA)
		maskMenu.AppendSeparator()
		tweak = wx.MenuItem(maskMenu, ID_MASK_TWEAK, 'Adjust Masking Parameters')
		maskMenu.AppendItem(tweak)
		
		## Bandpass Menu
		bandpassMenu.AppendRadioItem(ID_BANDPASS_OFF, 'Off')
		bandpassMenu.AppendRadioItem(ID_BANDPASS_ON,  'On')
		bandpassMenu.AppendSeparator()
		recompute = wx.MenuItem(bandpassMenu, ID_BANDPASS_RECOMPUTE, 'Recompute Fits')
		bandpassMenu.AppendItem(recompute)
		
		## Details Menu
		cf = wx.MenuItem(detailsMenu, ID_DETAIL_SUMMARY, 'Current File Info.')
		detailsMenu.AppendItem(cf)
		self.examineFileButton = wx.MenuItem(detailsMenu, ID_DETAIL_SUBFILE, 'Examine File')
		detailsMenu.AppendItem(self.examineFileButton)
		detailsMenu.AppendSeparator()
		zm = wx.MenuItem(detailsMenu, ID_DETAIL_WATERFALL, 'Zoomable Waterfall')
		detailsMenu.AppendItem(zm)
		zd = wx.MenuItem(detailsMenu, ID_DETAIL_DRIFTCURVE, 'Zoomable Drift Curve')
		detailsMenu.AppendItem(zd)

		# Creating the menubar.
		menuBar.Append(fileMenu,"&File") # Adding the "filemenu" to the MenuBar
		menuBar.Append(colorMenu, "&Color")
		menuBar.Append(dataMenu, "&Data")
		menuBar.Append(maskMenu, "&Mask")
		menuBar.Append(bandpassMenu, "&Bandpass")
		menuBar.Append(detailsMenu, "D&etails")
		self.SetMenuBar(menuBar)  # Adding the MenuBar to the Frame content.
		
		vbox = wx.BoxSizer(wx.VERTICAL)
		
		# Add waterfall plot
		panel1 = wx.Panel(self, -1)
		hbox1 = wx.BoxSizer(wx.HORIZONTAL)
		self.figure1a = Figure(figsize=(2,2))
		self.canvas1a = FigureCanvasWxAgg(panel1, -1, self.figure1a)
		hbox1.Add(self.canvas1a, 3, wx.ALIGN_LEFT | wx.EXPAND)
		
		# Add a saturation fraction plot
		self.figure1b = Figure()
		self.canvas1b = FigureCanvasWxAgg(panel1, -1, self.figure1b)
		hbox1.Add(self.canvas1b, 1, wx.ALIGN_CENTER | wx.EXPAND)
		
		# Add a total power with time plot
		self.figure1c = Figure(figsize=(2,2))
		self.canvas1c = FigureCanvasWxAgg(panel1, -1, self.figure1c)
		hbox1.Add(self.canvas1c, 3, wx.ALIGN_RIGHT | wx.EXPAND)
		panel1.SetSizer(hbox1)
		vbox.Add(panel1, 1, wx.EXPAND)
		
		# Add a spectrum plot (with toolbar)
		panel3 = wx.Panel(self, -1)
		hbox3 = wx.BoxSizer(wx.VERTICAL)
		self.figure2 = Figure()
		self.canvas2 = FigureCanvasWxAgg(panel3, -1, self.figure2)
		self.toolbar = NavigationToolbar2WxAgg(self.canvas2)
		self.toolbar.Realize()
		hbox3.Add(self.canvas2, 1, wx.EXPAND)
		hbox3.Add(self.toolbar, 0, wx.ALIGN_LEFT | wx.FIXED_MINSIZE)
		panel3.SetSizer(hbox3)
		vbox.Add(panel3, 1, wx.EXPAND)
		
		# Use some sizers to see layout options
		self.SetSizer(vbox)
		self.SetAutoLayout(1)
		vbox.Fit(self)
		
	def initEvents(self):
		self.Bind(wx.EVT_MENU, self.onOpen, id=ID_OPEN)
		self.Bind(wx.EVT_MENU, self.onSave, id=ID_SAVE)
		self.Bind(wx.EVT_MENU, self.onSaveAs, id=ID_SAVE_AS)
		self.Bind(wx.EVT_MENU, self.onExit, id=ID_QUIT)
		
		self.Bind(wx.EVT_MENU, self.onAutoscale, id=ID_COLOR_AUTO)
		self.Bind(wx.EVT_MENU, self.onAdjust, id=ID_COLOR_ADJUST)
		
		self.Bind(wx.EVT_MENU, self.onTuning1product1, id=ID_TUNING1_1)
		self.Bind(wx.EVT_MENU, self.onTuning1product2, id=ID_TUNING1_2)
		self.Bind(wx.EVT_MENU, self.onTuning1product3, id=ID_TUNING1_3)
		self.Bind(wx.EVT_MENU, self.onTuning1product4, id=ID_TUNING1_4)
		self.Bind(wx.EVT_MENU, self.onTuning2product1, id=ID_TUNING2_1)
		self.Bind(wx.EVT_MENU, self.onTuning2product2, id=ID_TUNING2_2)
		self.Bind(wx.EVT_MENU, self.onTuning2product3, id=ID_TUNING2_3)
		self.Bind(wx.EVT_MENU, self.onTuning2product4, id=ID_TUNING2_4)
		self.Bind(wx.EVT_MENU, self.onRangeChange, id=ID_CHANGE_RANGE)
		self.Bind(wx.EVT_MENU, self.onObservationChange, id=ID_CHANGE_OBSID)
		
		self.Bind(wx.EVT_MENU, self.onMaskSuggestCurrent, id=ID_MASK_SUGGEST_CURRENT)
		self.Bind(wx.EVT_MENU, self.onMaskSuggestAll, id=ID_MASK_SUGGEST_ALL)
		self.Bind(wx.EVT_MENU, self.onMaskResetCurrent, id=ID_MASK_RESET_CURRENT)
		self.Bind(wx.EVT_MENU, self.onMaskResetAll, id=ID_MASK_RESET_ALL)
		self.Bind(wx.EVT_MENU, self.onMaskTweak, id=ID_MASK_TWEAK)
		
		self.Bind(wx.EVT_MENU, self.onBandpassOn, id=ID_BANDPASS_ON)
		self.Bind(wx.EVT_MENU, self.onBandpassOff, id=ID_BANDPASS_OFF)
		self.Bind(wx.EVT_MENU, self.onBandpassRecompute, id=ID_BANDPASS_RECOMPUTE)
		
		self.Bind(wx.EVT_MENU, self.onFileDetails, id=ID_DETAIL_SUMMARY)
		self.Bind(wx.EVT_MENU, self.onExamineFile, id=ID_DETAIL_SUBFILE)
		self.Bind(wx.EVT_MENU, self.onZoomWaterfall, id=ID_DETAIL_WATERFALL)
		self.Bind(wx.EVT_MENU, self.onZoomDrift, id=ID_DETAIL_DRIFTCURVE)
		
		# Key events
		self.canvas1a.Bind(wx.EVT_KEY_UP, self.onKeyPress)
		self.canvas1b.Bind(wx.EVT_KEY_UP, self.onKeyPress)
		self.canvas2.Bind(wx.EVT_KEY_UP,  self.onKeyPress)
		
		# Make the plots resizable
		self.Bind(wx.EVT_PAINT, self.resizePlots)
		self.Bind(wx.EVT_SIZE, self.onSize)
	
	def onOpen(self, event):
		"""
		Open a file.
		"""
		
		dlg = wx.FileDialog(self, "Choose a file", self.dirname, "", "HDF5 (*.hdf5;*.h5)|*.hdf5;*.h5|All Files|*.*", wx.OPEN)
		if dlg.ShowModal() == wx.ID_OK:
			self.filename = dlg.GetFilename()
			self.dirname = dlg.GetDirectory()
			self.data = Waterfall_GUI(self)
			
			# Display the time range dialog box and update the image
			self.onRangeChange(None)
			
			if self.cAdjust is not None:
				try:
					self.cAdjust.Close()
				except:
					pass
				self.cAdjust = None
		dlg.Destroy()
		
		if self.data.filename is None:
			self.examineFileButton.Enable(False)
		else:
			self.examineFileButton.Enable(True)
		
	def onSave(self, event):
		"""
		Save the data mask to a new NPZ file.
		"""
		
		if self.data.filename == '':
			self.onSaveAs(event)
		else:
			h = h5py.File(self.data.filename, 'a')
			obs = h.get('Observation%i' % self.data.obsID, None)
			
			tuning1 = obs.get('Tuning1', None)
			tuning2 = obs.get('Tuning2', None)
			
			dataProducts = list(tuning1)
			mapper = {'XX': 0, 'I': 0, 'XY': 1, 'Q': 1, 'YX': 2, 'U': 2, 'YY': 3, 'V': 3}
			for exclude in ('freq', 'Saturation', 'SpectralKurtosis'):
				try:
					ind = dataProducts.index(exclude)
					del dataProducts[ind]
				except ValueError:
					pass
			
			o = self.data.iOffset
			d = self.data.iDuration
			
			mask1 = tuning1.get('Mask', None)
			mask2 = tuning2.get('Mask', None)
			
			if mask1 is None:
				mask1 = tuning1.create_group('Mask')
				for p in dataProducts:
					mask1.create_dataset(p, tuning1[p].shape, 'bool')
			for p in dataProducts:
				ind = 4*(1-1) + mapper[p]
				mask1P = mask1.get(p, None)
				mask1P[o:o+d,:] = self.data.spec.mask[:,ind,:]
				
			if mask2 is None:
				mask2 = tuning2.create_group('Mask')
				for p in dataProducts:
					mask2.create_dataset(p, tuning1[p].shape, 'bool')
			for p in dataProducts:
				ind = 4*(2-1) + mapper[p]
				mask2P = mask2.get(p, None)
				mask2P[o:o+d,:] = self.data.spec.mask[:,ind,:]
			
			h.close()

	def onSaveAs(self, event):
		"""
		Save the current observation to a new SD file.
		"""
		
		dialog = wx.FileDialog(self, "Select Output File", self.dirname, '', 'HDF5 Files (*.hdf5)|*.hdf5|All Files (*.*)|*.*', wx.SAVE|wx.FD_OVERWRITE_PROMPT)
			
		if dialog.ShowModal() == wx.ID_OK:
			self.dirname = dialog.GetDirectory()
			self.filename = dialog.GetPath()
			
			hOld = h5py.File(self.data.filename, 'r')
			hNew = h5py.File(self.filename, 'w')
			
			for name in hOld.attrs.keys():
				hNew.attrs[name] = hOld.attrs[name]
			
			for name in hOld.keys():
				hOld.copy(name, hNew)
				
			hOld.close()
			
			obs = hNew.get('Observation%i' % self.data.obsID, None)
			
			tuning1 = obs.get('Tuning1', None)
			tuning2 = obs.get('Tuning2', None)
			
			dataProducts = list(tuning1)
			mapper = {'XX': 0, 'I': 0, 'XY': 1, 'Q': 1, 'YX': 2, 'U': 2, 'YY': 3, 'V': 3}
			for exclude in ('freq', 'Saturation', 'SpectralKurtosis'):
				try:
					ind = dataProducts.index(exclude)
					del dataProducts[ind]
				except ValueError:
					pass
			
			o = self.data.iOffset
			d = self.data.iDuration
			
			mask1 = tuning1.get('Mask', None)
			mask2 = tuning2.get('Mask', None)
			
			if mask1 is None:
				mask1 = tuning1.create_group('Mask')
				for p in dataProducts:
					mask1.create_dataset(p, tuning1[p].shape, 'bool')
			for p in dataProducts:
				ind = 4*(1-1) + mapper[p]
				mask1P = mask1.get(p, None)
				mask1P[o:o+d,:] = self.data.spec.mask[:,ind,:]
				
			if mask2 is None:
				mask2 = tuning2.create_group('Mask')
				for p in dataProducts:
					mask2.create_dataset(p, tuning1[p].shape, 'bool')
			for p in dataProducts:
				ind = 4*(2-1) + mapper[p]
				mask2P = mask2.get(p, None)
				mask2P[o:o+d,:] = self.data.spec.mask[:,ind,:]
			
			hNew.close()
			
		dialog.Destroy()
		
	def onExit(self, event):
		"""
		Quit plotWaterfall.
		"""
		
		self.Close(True)
		
	def onAutoscale(self, event):
		"""
		Auto-scale the current data display.
		"""
		
		wx.BeginBusyCursor()
		wx.Yield()
		
		i = self.data.index
		toUse = numpy.arange(self.data.spec.shape[2]/10, 9*self.data.spec.shape[2]/10)
		
		try:
			from _helper import FastAxis1Percentiles5And99
			if self.data.bandpass:
				self.data.limitsBandpass[i] = FastAxis1Percentiles5And99(self.data.specBandpass.data, i, chanMin=self.data.spec.shape[2]/10, chanMax=9*self.data.spec.shape[2]/10)
			else:
				self.data.limits[i] = FastAxis1Percentiles5And99(self.data.spec.data, i)
		except ImportError:
			if self.data.bandpass:
				self.data.limitsBandpass[i] = [percentile(self.data.specBandpass[:,i,toUse].ravel(), 5), percentile(self.data.specBandpass[:,i,toUse].ravel(), 99)] 
			else:
				self.data.limits[i] = [percentile(self.data.spec[:,i,:].ravel(), 5), percentile(self.data.spec[:,i,:].ravel(), 99)]
			
		if self.data.usedB:
			if self.data.bandpass:
				self.data.limitsBandpass[i] = [to_dB(v) for v in self.data.limitsBandpass[i]]
			else:
				self.data.limits[i] = [to_dB(v) for v in self.data.limits[i]]
				
		self.data.draw()
		self.data.drawSpectrum(self.data.spectrumClick)
		self.data.makeMark(self.data.spectrumClick)
		
		wx.EndBusyCursor()
		
	def onAdjust(self, event):
		"""
		Bring up the colorbar adjustment dialog window.
		"""
		
		ContrastAdjust(self)
		
	def onTuning1product1(self, event):
		"""
		Display tuning 1, data product 1 (XX or I)
		"""
		
		wx.BeginBusyCursor()
		wx.Yield()
		
		self.data.index = 0
		self.data.draw()
		if self.data.spectrumClick is not None:
			self.data.drawSpectrum(self.data.spectrumClick)
			
		wx.EndBusyCursor()
		
	def onTuning1product2(self, event):
		"""
		Display tuning 1, data product 1 (XY or Q)
		"""
		
		wx.BeginBusyCursor()
		wx.Yield()
		
		self.data.index = 1
		self.data.draw()
		if self.data.spectrumClick is not None:
			self.data.drawSpectrum(self.data.spectrumClick)
			
		wx.EndBusyCursor()
		
	def onTuning1product3(self, event):
		"""
		Display tuning 1, data product 3 (YX or U)
		"""
		
		wx.BeginBusyCursor()
		wx.Yield()
		
		self.data.index = 2
		self.data.draw()
		if self.data.spectrumClick is not None:
			self.data.drawSpectrum(self.data.spectrumClick)
			
		wx.EndBusyCursor()
		
	def onTuning1product4(self, event):
		"""
		Display tuning 1, data product 4 (YY or V)
		"""
		
		wx.BeginBusyCursor()
		wx.Yield()
		
		self.data.index = 3
		self.data.draw()
		if self.data.spectrumClick is not None:
			self.data.drawSpectrum(self.data.spectrumClick)
			
		wx.EndBusyCursor()
		
	def onTuning2product1(self, event):
		"""
		Display tuning 2, data product 1 (XX or I)
		"""
		
		wx.BeginBusyCursor()
		wx.Yield()
		
		self.data.index = 4
		self.data.draw()
		if self.data.spectrumClick is not None:
			self.data.drawSpectrum(self.data.spectrumClick)
			
		wx.EndBusyCursor()
		
	def onTuning2product2(self, event):
		"""
		Display tuning 2, data product 1 (XY or Q)
		"""
		
		wx.BeginBusyCursor()
		wx.Yield()
		
		self.data.index = 5
		self.data.draw()
		if self.data.spectrumClick is not None:
			self.data.drawSpectrum(self.data.spectrumClick)
			
		wx.EndBusyCursor()
		
	def onTuning2product3(self, event):
		"""
		Display tuning 2, data product 3 (YX or U)
		"""
		
		wx.BeginBusyCursor()
		wx.Yield()
		
		self.data.index = 6
		self.data.draw()
		if self.data.spectrumClick is not None:
			self.data.drawSpectrum(self.data.spectrumClick)
			
		wx.EndBusyCursor()
		
	def onTuning2product4(self, event):
		"""
		Display tuning 2, data product 4 (YY or V)
		"""
		
		wx.BeginBusyCursor()
		wx.Yield()
		
		self.data.index = 7
		self.data.draw()
		if self.data.spectrumClick is not None:
			self.data.drawSpectrum(self.data.spectrumClick)
			
		wx.EndBusyCursor()
		
	def onRangeChange(self, event):
		"""
		Display a dialog box to change the time range displayed.
		"""
		
		if event is None:
			mode = 'New'
		else:
			mode = 'Adjust'
			
		TimeRangeAdjust(self, mode=mode)
		
	def onObservationChange(self, event):
		"""
		Display a dialog box to change the observation currently being 
		displayed.
		"""
		
		SwitchObservation(self)
		
	def onMaskSuggestCurrent(self, event):
		"""
		Suggest a series of frequency and time-based masks to apply to 
		the current tuning/polarization.
		"""
		
		wx.BeginBusyCursor()
		wx.Yield()
		
		self.data.suggestMask(self.data.index)
			
		self.data.draw()
		self.data.drawSpectrum(self.data.spectrumClick)
		self.data.makeMark(self.data.spectrumClick)
		
		wx.EndBusyCursor()
	
	def onMaskSuggestAll(self, event):
		"""
		Suggest a series of frequency and time-based masks to apply to
		all data streams.
		"""
		
		wx.BeginBusyCursor()
		wx.Yield()
		
		for i in xrange(self.data.spec.shape[1]):
			self.data.suggestMask(i)
			
		self.data.draw()
		self.data.drawSpectrum(self.data.spectrumClick)
		self.data.makeMark(self.data.spectrumClick)
		
		wx.EndBusyCursor()
	
	def onMaskResetCurrent(self, event):
		"""
		Reset the current mask.
		"""
		
		wx.BeginBusyCursor()
		wx.Yield()
		
		self.data.resetMask(self.data.index)
		
		self.data.draw()
		self.data.drawSpectrum(self.data.spectrumClick)
		self.data.makeMark(self.data.spectrumClick)
		
		wx.EndBusyCursor()
		
	def onMaskResetAll(self, event):
		"""
		Reset all masks.
		"""
		
		wx.BeginBusyCursor()
		wx.Yield()
		
		for i in xrange(self.data.spec.shape[1]):
			self.data.resetMask(i)
		
		self.data.draw()
		self.data.drawSpectrum(self.data.spectrumClick)
		self.data.makeMark(self.data.spectrumClick)
		
		wx.EndBusyCursor()
		
	def onMaskTweak(self, event):
		"""
		Tweak the masking parameters.
		"""
		
		MaskingAdjust(self)
	
	def onBandpassOn(self, event):
		"""
		Enable bandpass correction.
		"""
		
		wx.BeginBusyCursor()
		wx.Yield()
		
		self.data.bandpass = True
		
		self.data.draw()
		self.data.drawSpectrum(self.data.spectrumClick)
		self.data.makeMark(self.data.spectrumClick)
		
		wx.EndBusyCursor()
		
	def onBandpassOff(self, event):
		"""
		Disable bandpass correction.
		"""
		
		wx.BeginBusyCursor()
		wx.Yield()
		
		self.data.bandpass = False
		
		self.data.draw()
		self.data.drawSpectrum(self.data.spectrumClick)
		self.data.makeMark(self.data.spectrumClick)
		
		wx.EndBusyCursor()
		
	def onBandpassRecompute(self, event):
		"""
		Recompute the bandpass fits and redraw the plots if needed.
		"""
		
		wx.BeginBusyCursor()
		wx.Yield()
		
		self.data.computeBandpass()
		
		if self.data.bandpass:
			self.data.draw()
			self.data.drawSpectrum(self.data.spectrumClick)
			self.data.makeMark(self.data.spectrumClick)
			
		wx.EndBusyCursor()
		
	def onFileDetails(self, event):
		"""
		Display a small info. window about the current file.
		"""
		
		# Get some basic parameter
		filename = self.data.filename
		beam = self.data.beam
		srate, sunit = bestFreqUnits(self.data.srate)
		tInt = self.data.tInt
		nInt = self.data.spec.shape[0]
		isAggregate = False if self.data.filenames is None else True
		tIntOrg = self.data.tIntOriginal
		tIntAct = self.data.tIntActual
		
		target = self.data.target
		if self.data.raUnit.lower().find('h') != -1:
			ra = ephem.hours(self.data.ra*numpy.pi/12.0)
		else:
			ra = ephem.hours(self.data.ra*numpy.pi/180.0)
		if self.data.decUnit.lower().find('h') != -1:
			dec = ephem.degrees(self.data.dec*numpy.pi/12.0)
		else:
			dec = ephem.degrees(self.data.dec*numpy.pi/180.0)
		mode = self.data.mode
		
		rbw = self.data.rbw
		rbwu = self.data.rbwUnit
		
		# Build the message string
		outString = """Filename: %s

Target Name: %s
RA: %s
Dec: %s
Observing Mode: %s

Beam:  %i
Sample Rate: %.3f %s

Integration Time:  %.3f seconds
RBW: %.3f %s
Number of Integrations:  %i

Aggregate File?  %s""" % (filename, target, ra, dec, mode, beam, srate, sunit, tInt, rbw, rbwu, nInt, isAggregate)

		# Expound on aggregate files
		if isAggregate:
			outString = """%s
Number of files contained:  %i
Original Integration Time:  %.3f seconds
Actual Integration Time:  %.3f seconds""" % (outString, len(self.data.filenames), tIntOrg, tIntAct)
		
		# Show the box
		box = wx.MessageDialog(self, outString, "File Details")
		box.ShowModal()
		
	def onExamineFile(self, event):
		"""
		For aggregated data sets, open a new plotWaterfall to look at the
		current file.
		"""
		
		# Make sure there is not another window running already.
		if self.examineWindow is not None:
			self.examineWindow.poll()
			if self.examineWindow.returncode is None:
				print "ERROR: another sub-file examination window is already running"
				return False
			else:
				self.examineWindow = None
		
		# Make sure we actually have a aggregated file first
		if self.data.filenames is None:
			print "ERROR: current file is not an aggregated file"
			return False
		
		# Make sure we have clicked
		try:
			dataY = int(round(self.data.spectrumClick / self.data.tInt))
		except:
			print "ERROR: no sub-file currently selected"
			return False
			
		# Make sure the target file exists
		filename = self.data.filenames[dataY]
		if not os.path.exists(filename):
			print "ERROR: cannot find file '%s', trying NPZ path" % filename
			basepath, junk = os.path.split(self.data.filename)
			junk, filename = os.path.split(self.data.filenames[dataY])
			filename = os.path.join(basepath, filename)
			
			if not os.path.exists(filename):
				print "ERROR: cannot find file '%s'" % filename
				return False
		
		# Get the current plotWaterfall.py path
		script = sys.argv[0]
		
		# Go
		self.examineWindow = subprocess.Popen([sys.executable, script, filename])
		
	def onZoomWaterfall(self, event):
		"""
		Create a zoomable waterfall plot window.
		"""
		
		WaterfallDisplay(self)
	
	def onZoomDrift(self, event):
		"""
		Create a zoomable drift curve plot window.
		"""
		
		DriftCurveDisplay(self)
		
	def onKeyPress(self, event):
		"""
		Move the current spectra mark up and down with a keypress.
		"""
		
		keycode = event.GetKeyCode()
		
		if keycode == wx.WXK_UP:
			## Move forward one integration time
			if self.data.spectrumClick is not None:
				self.data.drawSpectrum(self.data.spectrumClick + self.data.tInt)
				self.data.makeMark(self.data.spectrumClick + self.data.tInt)
		elif keycode == wx.WXK_DOWN:
			## Move backward one integration time
			if self.data.spectrumClick is not None:
				self.data.drawSpectrum(self.data.spectrumClick - self.data.tInt)
				self.data.makeMark(self.data.spectrumClick - self.data.tInt)
		elif keycode == wx.WXK_SPACE:
			## Mask the current integration
			if self.data.spectrumClick is not None:
				dataY = int(round(self.data.spectrumClick / self.data.tInt ))
				
				self.data.spec.mask[dataY, self.data.index, :] = True
				self.data.specBandpass.mask[dataY, self.data.index, :] = True
				self.data.timeMask[dataY, self.data.index] = True
				
				self.data.draw()
		elif keycode == 85:
			## Unmask the current integration
			if self.data.spectrumClick is not None:
				dataY = int(round(self.data.spectrumClick / self.data.tInt ))
				
				self.data.spec.mask[dataY, self.data.index, :] = self.data.freqMask[self.data.index,:]
				self.data.specBandpass.mask[dataY, self.data.index, :] = self.data.freqMask[self.data.index,:]
				self.data.timeMask[dataY, self.data.index] = False
				
				self.data.draw()
		else:
			pass
			
		event.Skip()
		
	def onSize(self, event):
		event.Skip()
		wx.CallAfter(self.resizePlots)
		
	def resizePlots(self, event=None):
		# Get the current size of the window and the navigation toolbar
		w, h = self.GetClientSize()
		wt, ht = self.toolbar.GetSize()
		
		# Come up with new figure sizes in inches.
		# NOTE:  The height of the spectrum plot at the bottom needs to be
		#        corrected for the height of the toolbar
		dpi = self.figure1a.get_dpi()
		newW0 = 1.0*w/dpi
		newW1 = 1.0*(3.*w/7)/dpi
		newW2 = 1.0*(1.*w/7)/dpi
		newW3 = 1.0*(3.*w/7)/dpi
		newH0 = 1.0*(h/2)/dpi
		newH1 = 1.0*(h/2-ht)/dpi
		
		# Set the figure sizes and redraw
		self.figure1a.set_size_inches((newW1, newH0))
		try:
			self.figure1a.tight_layout()
		except:
			pass
		self.figure1a.canvas.draw()
		
		self.figure1b.set_size_inches((newW2, newH0))
		try:
			self.figure1b.tight_layout()
		except:
			pass
		self.figure1b.canvas.draw()
		
		self.figure1c.set_size_inches((newW3, newH0))
		try:
			self.figure1c.tight_layout()
		except:
			pass
		self.figure1c.canvas.draw()
		
		self.figure2.set_size_inches((newW0, newH1))
		try:
			self.figure2.tight_layout()
		except:
			pass
		self.figure2.canvas.draw()
		
	def GetToolBar(self):
		# You will need to override GetToolBar if you are using an 
		# unmanaged toolbar in your frame
		return self.toolbar


ID_CONTRAST_UPR_INC = 100
ID_CONTRAST_UPR_DEC = 101
ID_CONTRAST_LWR_INC = 102
ID_CONTRAST_LWR_DEC = 103
ID_CONTRAST_OK = 104

class ContrastAdjust(wx.Frame):
	def __init__ (self, parent):	
		wx.Frame.__init__(self, parent, title='Contrast Adjustment', size=(330, 175))
		
		self.parent = parent
		
		self.initUI()
		self.initEvents()
		self.Show()
		
		self.parent.cAdjust = self
		
	def initUI(self):
		row = 0
		panel = wx.Panel(self)
		sizer = wx.GridBagSizer(5, 5)
		
		if self.parent.data.bandpass:
			typ = wx.StaticText(panel, label='Tuning %i, Pol. %i - Bandpass' % (self.parent.data.index/2+1, self.parent.data.index%2))
		else:
			typ = wx.StaticText(panel, label='Tuning %i, Pol. %i' % (self.parent.data.index/2+1, self.parent.data.index%2))
		sizer.Add(typ, pos=(row+0, 0), span=(1,4), flag=wx.EXPAND|wx.LEFT|wx.RIGHT, border=5)
		
		row += 1
		
		upr = wx.StaticText(panel, label='Upper Limit:')
		uprText = wx.TextCtrl(panel)
		uprDec = wx.Button(panel, ID_CONTRAST_UPR_DEC, '-', size=(56, 28))
		uprInc = wx.Button(panel, ID_CONTRAST_UPR_INC, '+', size=(56, 28))
		
		lwr = wx.StaticText(panel, label='Lower Limit:')
		lwrText = wx.TextCtrl(panel)
		lwrDec = wx.Button(panel, ID_CONTRAST_LWR_DEC, '-', size=(56, 28))
		lwrInc = wx.Button(panel, ID_CONTRAST_LWR_INC, '+', size=(56, 28))
		
		rng = wx.StaticText(panel, label='Range:')
		rngText = wx.TextCtrl(panel, style=wx.TE_READONLY)
		
		sizer.Add(upr,     pos=(row+0, 0), span=(1, 1), flag=wx.EXPAND|wx.LEFT|wx.RIGHT, border=5)
		sizer.Add(uprText, pos=(row+0, 1), span=(1, 1), flag=wx.EXPAND|wx.LEFT|wx.RIGHT, border=5)
		sizer.Add(uprDec,  pos=(row+0, 2), span=(1, 1), flag=wx.EXPAND|wx.LEFT|wx.RIGHT, border=5)
		sizer.Add(uprInc,  pos=(row+0, 3), span=(1, 1), flag=wx.EXPAND|wx.LEFT|wx.RIGHT, border=5)
		sizer.Add(lwr,     pos=(row+1, 0), span=(1, 1), flag=wx.EXPAND|wx.LEFT|wx.RIGHT, border=5)
		sizer.Add(lwrText, pos=(row+1, 1), span=(1, 1), flag=wx.EXPAND|wx.LEFT|wx.RIGHT, border=5)
		sizer.Add(lwrDec,  pos=(row+1, 2), span=(1, 1), flag=wx.EXPAND|wx.LEFT|wx.RIGHT, border=5)
		sizer.Add(lwrInc,  pos=(row+1, 3), span=(1, 1), flag=wx.EXPAND|wx.LEFT|wx.RIGHT, border=5)
		sizer.Add(rng,     pos=(row+2, 0), span=(1, 1), flag=wx.EXPAND|wx.LEFT|wx.RIGHT, border=5)
		sizer.Add(rngText, pos=(row+2, 1), span=(1, 1), flag=wx.EXPAND|wx.LEFT|wx.RIGHT, border=5)
		
		line = wx.StaticLine(panel)
		sizer.Add(line, pos=(row+3, 0), span=(1, 4), flag=wx.EXPAND|wx.BOTTOM, border=10)
		
		row += 4
		
		#
		# Buttons
		#
		
		ok = wx.Button(panel, ID_CONTRAST_OK, 'Ok', size=(56, 28))
		sizer.Add(ok, pos=(row+0, 3), flag=wx.RIGHT|wx.BOTTOM, border=5)
		
		panel.SetSizerAndFit(sizer)

		self.uText = uprText
		self.lText = lwrText
		self.rText = rngText
		
		#
		# Set current values
		#
		index = self.parent.data.index
		if self.parent.data.bandpass:
			self.uText.SetValue('%.1f' % self.parent.data.limitsBandpass[index][1])
			self.lText.SetValue('%.1f' % self.parent.data.limitsBandpass[index][0])
			self.rText.SetValue('%.1f' % self.__getRange(index))
		else:
			self.uText.SetValue('%.1f' % self.parent.data.limits[index][1])
			self.lText.SetValue('%.1f' % self.parent.data.limits[index][0])
			self.rText.SetValue('%.1f' % self.__getRange(index))
		
	def initEvents(self):
		self.Bind(wx.EVT_BUTTON, self.onUpperDecrease, id=ID_CONTRAST_UPR_DEC)
		self.Bind(wx.EVT_BUTTON, self.onUpperIncrease, id=ID_CONTRAST_UPR_INC)
		self.Bind(wx.EVT_BUTTON, self.onLowerDecrease, id=ID_CONTRAST_LWR_DEC)
		self.Bind(wx.EVT_BUTTON, self.onLowerIncrease, id=ID_CONTRAST_LWR_INC)
		
		self.Bind(wx.EVT_BUTTON, self.onOk, id=ID_CONTRAST_OK)
		
		self.Bind(wx.EVT_KEY_UP, self.onKeyPress)
		
	def onKeyPress(self, event):
		keycode = event.GetKeyCode()
		
		if keycode == wx.WXK_RETURN:
			index = self.parent.data.index
			if self.parent.data.bandpass:
				self.parent.data.limitsBandpass[index][0] = float(self.lText.GetValue())
				self.parent.data.limitsBandpass[index][1] = float(self.uText.GetValue())
			else:
				self.parent.data.limits[index][0] = float(self.lText.GetValue())
				self.parent.data.limits[index][1] = float(self.uText.GetValue())
			self.rText.SetValue('%.1f' % self.__getRange(index))
			self.parent.data.draw()
		else:
			pass
		
	def onUpperDecrease(self, event):
		index = self.parent.data.index
		if self.parent.data.bandpass:
			self.parent.data.limitsBandpass[index][1] -= self.__getIncrement(index)
			self.uText.SetValue('%.1f' % self.parent.data.limitsBandpass[index][1])
		else:
			self.parent.data.limits[index][1] -= self.__getIncrement(index)
			self.uText.SetValue('%.1f' % self.parent.data.limits[index][1])
		self.rText.SetValue('%.1f' % self.__getRange(index))
		self.parent.data.draw()
		
	def onUpperIncrease(self, event):
		index = self.parent.data.index
		if self.parent.data.bandpass:
			self.parent.data.limitsBandpass[index][1] += self.__getIncrement(index)
			self.uText.SetValue('%.1f' % self.parent.data.limitsBandpass[index][1])
		else:
			self.parent.data.limits[index][1] += self.__getIncrement(index)
			self.uText.SetValue('%.1f' % self.parent.data.limits[index][1])
		self.rText.SetValue('%.1f' % self.__getRange(index))
		self.parent.data.draw()
		
	def onLowerDecrease(self, event):
		index = self.parent.data.index
		if self.parent.data.bandpass:
			self.parent.data.limitsBandpass[index][0] -= self.__getIncrement(index)
			self.lText.SetValue('%.1f' % self.parent.data.limitsBandpass[index][0])
		else:
			self.parent.data.limits[index][0] -= self.__getIncrement(index)
			self.lText.SetValue('%.1f' % self.parent.data.limits[index][0])
		self.rText.SetValue('%.1f' % self.__getRange(index))
		self.parent.data.draw()
		
	def onLowerIncrease(self, event):
		index = self.parent.data.index
		if self.parent.data.bandpass:
			self.parent.data.limitsBandpass[index][0] += self.__getIncrement(index)
			self.lText.SetValue('%.1f' % self.parent.data.limitsBandpass[index][0])
		else:
			self.parent.data.limits[index][0] += self.__getIncrement(index)
			self.lText.SetValue('%.1f' % self.parent.data.limits[index][0])
		self.rText.SetValue('%.1f' % self.__getRange(index))
		self.parent.data.draw()
		
	def onOk(self, event):
		self.parent.cAdjust = None
		self.Close()
	
	def __getRange(self, index):
		if self.parent.data.bandpass:
			return (self.parent.data.limitsBandpass[index][1] - self.parent.data.limitsBandpass[index][0])
		else:
			return (self.parent.data.limits[index][1] - self.parent.data.limits[index][0])
		
	def __getIncrement(self, index):
		return 0.1*self.__getRange(index)


ID_RANGE_WHOLE = 100
ID_RANGE_OK = 101
ID_RANGE_CANCEL = 102

class TimeRangeAdjust(wx.Frame):
	def __init__ (self, parent, mode='Adjust'):	
		wx.Frame.__init__(self, parent, title='Time Range to Display', size=(330, 175), style=wx.STAY_ON_TOP|wx.FRAME_FLOAT_ON_PARENT)
		
		self.parent = parent
		self.mode = mode
		
		self.initUI()
		self.initEvents()
		self.Show()
		
	def initUI(self):
		row = 0
		panel = wx.Panel(self)
		sizer = wx.GridBagSizer(5, 5)
		
		offL = wx.StaticText(panel, label='Offset from file start:')
		sizer.Add(offL, pos=(row+0, 0), span=(1, 2), flag=wx.EXPAND|wx.LEFT|wx.RIGHT, border=5)
		row += 1
		
		self.offsetText = wx.TextCtrl(panel)
		offU = wx.StaticText(panel, label='seconds')
		sizer.Add(self.offsetText, pos=(row+0, 0), span=(1, 1), flag=wx.EXPAND|wx.LEFT|wx.RIGHT, border=5)
		sizer.Add(offU, pos=(row+0, 1), span=(1, 1), flag=wx.EXPAND|wx.LEFT|wx.RIGHT, border=5)
		row += 1
		
		durL = wx.StaticText(panel, label='Duration after offset:')
		sizer.Add(durL, pos=(row+0, 0), span=(1, 2), flag=wx.EXPAND|wx.LEFT|wx.RIGHT, border=5)
		row += 1
		
		self.durationText = wx.TextCtrl(panel)
		durU = wx.StaticText(panel, label='seconds')
		sizer.Add(self.durationText, pos=(row+0, 0), span=(1, 1), flag=wx.EXPAND|wx.LEFT|wx.RIGHT, border=5)
		sizer.Add(durU, pos=(row+0, 1), span=(1, 1), flag=wx.EXPAND|wx.LEFT|wx.RIGHT, border=5)
		row += 1
		
		self.wholeFileButton = wx.CheckBox(panel, ID_RANGE_WHOLE, 'Display whole observation')
		sizer.Add(self.wholeFileButton, pos=(row+0, 0), span=(1, 2), flag=wx.EXPAND|wx.LEFT|wx.RIGHT, border=5)
		row += 1
		
		#
		# Buttons
		#
		
		ok = wx.Button(panel, ID_RANGE_OK, 'Ok', size=(56, 28))
		sizer.Add(ok, pos=(row+0, 0), flag=wx.RIGHT|wx.BOTTOM, border=5)
		cancel = wx.Button(panel, ID_RANGE_CANCEL, 'Cancel', size=(56, 28))
		sizer.Add(cancel, pos=(row+0, 1),  flag=wx.RIGHT|wx.BOTTOM, border=5)
		
		#
		# Fill in values
		#
		self.offsetText.SetValue("%.3f" % self.parent.offset)
		self.durationText.SetValue("%.3f" % self.parent.duration)
		if self.parent.offset == 0 and self.parent.duration < 0:
			self.wholeFileButton.SetValue(True)
			self.offsetText.Disable()
			self.durationText.Disable()
			
		#
		# Set the operational mode
		#
		if self.mode != 'Adjust':
			cancel.Disable()
		
		panel.SetSizerAndFit(sizer)
		
	def initEvents(self):
		self.Bind(wx.EVT_CHECKBOX, self.onWholeFileToggle, id=ID_RANGE_WHOLE)
		
		self.Bind(wx.EVT_BUTTON, self.onOk, id=ID_RANGE_OK)
		self.Bind(wx.EVT_BUTTON, self.onCancel, id=ID_RANGE_CANCEL)
		
	def onWholeFileToggle(self, event):
		if self.wholeFileButton.GetValue():
			self.offsetText.Disable()
			self.durationText.Disable()
		else:
			self.offsetText.Enable()
			self.durationText.Enable()
		
	def onOk(self, event):
		# Get values
		if self.wholeFileButton.GetValue():
			newOffset = 0.0
			newDuration = -1.0
		else:
			newOffset = float(self.offsetText.GetValue())
			newDuration = float(self.durationText.GetValue())
		
		# Figure out if we need to update
		needToUpdate = False
		if self.parent.offset != newOffset:
			needToUpdate = True
			self.parent.offset = newOffset
		if self.parent.duration != newDuration:
			needToUpdate = True
			self.parent.duration = newDuration
		
		try:
			if needToUpdate or self.mode == 'New':
				self.parent.data.loadData(os.path.join(self.parent.dirname, self.parent.filename))
				self.parent.data.render()
				self.parent.data.draw()
		except Exception, e:
			print "ERROR: %s" % str(e)
		else:
			self.Close()
		
	def onCancel(self, event):
		self.Close()


ID_OBSID_OK = 101
ID_OBSID_CANCEL = 102

class SwitchObservation(wx.Frame):
	def __init__ (self, parent):	
		wx.Frame.__init__(self, parent, title='Observation to Display', size=(330, 175), style=wx.STAY_ON_TOP|wx.FRAME_FLOAT_ON_PARENT)
		
		self.parent = parent
		
		self.initUI()
		self.initEvents()
		
		h,w = self.GetEffectiveMinSize()
		self.SetSize((h,w))
		self.Fit()
		
		self.Show()
		
	def initUI(self):
		obsList = self.getCurrentObsList()
		
		row = 0
		panel = wx.Panel(self)
		sizer = wx.GridBagSizer(5, 5)
		
		self.rbList = []
		for i,obs in enumerate(obsList):
			if i == 0:
				rb = wx.RadioButton(panel, -1, '%s in mode %s (#%i)' % (obs[1], obs[2], obs[0]), style=wx.RB_GROUP)
			else:
				rb = wx.RadioButton(panel, -1, '%s in mode %s (#%i)' % (obs[1], obs[2], obs[0]))
				
			sizer.Add(rb, pos=(row+0, 0), flag=wx.EXPAND|wx.LEFT|wx.RIGHT, border=5)
			row += 1
			self.rbList.append((rb,obs[0]))
			
		#
		# Buttons
		#
		
		ok = wx.Button(panel, ID_OBSID_OK, 'Ok', size=(56, 28))
		sizer.Add(ok, pos=(row+0, 0), flag=wx.RIGHT|wx.BOTTOM, border=5)
		cancel = wx.Button(panel, ID_OBSID_CANCEL, 'Cancel', size=(56, 28))
		sizer.Add(cancel, pos=(row+0, 1),  flag=wx.RIGHT|wx.BOTTOM, border=5)
		
		#
		# Fill in values
		#
		currObsID = self.parent.data.obsID
		for rb,id in self.rbList:
			if id == currObsID:
				rb.SetValue(True)
			else:
				rb.SetValue(False)
				
		panel.SetSizerAndFit(sizer)
		
	def initEvents(self):
		self.Bind(wx.EVT_BUTTON, self.onOk, id=ID_OBSID_OK)
		self.Bind(wx.EVT_BUTTON, self.onCancel, id=ID_OBSID_CANCEL)
		
	def onOk(self, event):
		# Get values
		for rb,id in self.rbList:
			if rb.GetValue():
				newID = id
				
		try:
			if newID != self.parent.data.obsID:
				self.parent.data.loadData(os.path.join(self.parent.dirname, self.parent.filename), obsID=newID)
				self.parent.data.render()
				self.parent.data.draw()
		except Exception, e:
			print "ERROR: %s" % str(e)
		else:
			self.Close()
		
	def onCancel(self, event):
		self.Close()
		
	def getCurrentObsList(self):
		"""
		Return a list of tuples, one for each observation in the current file.
		"""
		
		obsList = []
		
		h = h5py.File(self.parent.data.filename)
		for obsID in list(h):
			id = int(obsID[11:])
			obs = h.get(obsID, None)
			obsList.append( (id, obs.attrs['TargetName'], obs.attrs['TrackingMode']) )
		h.close()
		
		return obsList


ID_BANDPASS_CUT_INC = 100
ID_BANDPASS_CUT_DEC = 101
ID_DRIFT_POLY_INC = 102
ID_DRIFT_POLY_DEC = 103
ID_DRIFT_CUT_INC  = 104
ID_DRIFT_CUT_DEC  = 105
ID_SK_SEC_INC = 106
ID_SK_SEC_DEC = 107
ID_SK_CUT_INC = 108
ID_SK_CUT_DEC = 109
ID_MASKING_OK = 110

class MaskingAdjust(wx.Frame):
	def __init__ (self, parent):	
		wx.Frame.__init__(self, parent, title='Masking Adjustment', size=(330, 350))
		
		self.parent = parent
		
		self.initUI()
		self.initEvents()
		self.Show()
		
		self.parent.cAdjust = self
		
	def initUI(self):
		row = 0
		panel = wx.Panel(self)
		sizer = wx.GridBagSizer(5, 5)
		
		typ = wx.StaticText(panel, label='Masking Parameters')
		sizer.Add(typ, pos=(row+0, 0), span=(1,4), flag=wx.EXPAND|wx.LEFT|wx.RIGHT, border=5)
		
		row += 1
		
		bp = wx.StaticText(panel, label='Bandpass Rention:')
		bpR = wx.StaticText(panel, label='Inner:')
		bpRText = wx.TextCtrl(panel, style=wx.TE_READONLY)
		bpRDec = wx.Button(panel, ID_BANDPASS_CUT_DEC, '-', size=(56, 28))
		bpRInc = wx.Button(panel, ID_BANDPASS_CUT_INC, '+', size=(56, 28))
		
		dc = wx.StaticText(panel, label='Drift Curve:')
		dcP = wx.StaticText(panel, label='Fit order:')
		dcPText = wx.TextCtrl(panel, style=wx.TE_READONLY)
		dcPDec = wx.Button(panel, ID_DRIFT_POLY_DEC, '-', size=(56, 28))
		dcPInc = wx.Button(panel, ID_DRIFT_POLY_INC, '+', size=(56, 28))
		
		dcC = wx.StaticText(panel, label='Threshold:')
		dcCText = wx.TextCtrl(panel, style=wx.TE_READONLY)
		dcCDec = wx.Button(panel, ID_DRIFT_CUT_DEC, '-', size=(56, 28))
		dcCInc = wx.Button(panel, ID_DRIFT_CUT_INC, '+', size=(56, 28))
		
		sk = wx.StaticText(panel, label='Spectral Kurtosis:')
		skS = wx.StaticText(panel, label='Sections:')
		skSText = wx.TextCtrl(panel, style=wx.TE_READONLY)
		skSDec = wx.Button(panel, ID_SK_SEC_DEC, '-', size=(56, 28))
		skSInc = wx.Button(panel, ID_SK_SEC_INC, '+', size=(56, 28))
		
		skC = wx.StaticText(panel, label='Threshold:')
		skCText = wx.TextCtrl(panel, style=wx.TE_READONLY)
		skCDec = wx.Button(panel, ID_SK_CUT_DEC, '-', size=(56, 28))
		skCInc = wx.Button(panel, ID_SK_CUT_INC, '+', size=(56, 28))
		
		sizer.Add(bp,      pos=(row+0, 0), span=(1, 4), flag=wx.EXPAND|wx.LEFT|wx.RIGHT, border=5)
		sizer.Add(bpR,     pos=(row+1, 0), span=(1, 1), flag=wx.EXPAND|wx.LEFT|wx.RIGHT, border=5)
		sizer.Add(bpRText, pos=(row+1, 1), span=(1, 1), flag=wx.EXPAND|wx.LEFT|wx.RIGHT, border=5)
		sizer.Add(bpRDec,  pos=(row+1, 2), span=(1, 1), flag=wx.EXPAND|wx.LEFT|wx.RIGHT, border=5)
		sizer.Add(bpRInc,  pos=(row+1, 3), span=(1, 1), flag=wx.EXPAND|wx.LEFT|wx.RIGHT, border=5)
		
		sizer.Add(dc,      pos=(row+2, 0), span=(1, 4), flag=wx.EXPAND|wx.LEFT|wx.RIGHT, border=5)
		sizer.Add(dcP,     pos=(row+3, 0), span=(1, 1), flag=wx.EXPAND|wx.LEFT|wx.RIGHT, border=5)
		sizer.Add(dcPText, pos=(row+3, 1), span=(1, 1), flag=wx.EXPAND|wx.LEFT|wx.RIGHT, border=5)
		sizer.Add(dcPDec,  pos=(row+3, 2), span=(1, 1), flag=wx.EXPAND|wx.LEFT|wx.RIGHT, border=5)
		sizer.Add(dcPInc,  pos=(row+3, 3), span=(1, 1), flag=wx.EXPAND|wx.LEFT|wx.RIGHT, border=5)
		sizer.Add(dcC,     pos=(row+4, 0), span=(1, 1), flag=wx.EXPAND|wx.LEFT|wx.RIGHT, border=5)
		sizer.Add(dcCText, pos=(row+4, 1), span=(1, 1), flag=wx.EXPAND|wx.LEFT|wx.RIGHT, border=5)
		sizer.Add(dcCDec,  pos=(row+4, 2), span=(1, 1), flag=wx.EXPAND|wx.LEFT|wx.RIGHT, border=5)
		sizer.Add(dcCInc,  pos=(row+4, 3), span=(1, 1), flag=wx.EXPAND|wx.LEFT|wx.RIGHT, border=5)
		
		sizer.Add(sk,      pos=(row+5, 0), span=(1, 4), flag=wx.EXPAND|wx.LEFT|wx.RIGHT, border=5)
		sizer.Add(skS,     pos=(row+6, 0), span=(1, 1), flag=wx.EXPAND|wx.LEFT|wx.RIGHT, border=5)
		sizer.Add(skSText, pos=(row+6, 1), span=(1, 1), flag=wx.EXPAND|wx.LEFT|wx.RIGHT, border=5)
		sizer.Add(skSDec,  pos=(row+6, 2), span=(1, 1), flag=wx.EXPAND|wx.LEFT|wx.RIGHT, border=5)
		sizer.Add(skSInc,  pos=(row+6, 3), span=(1, 1), flag=wx.EXPAND|wx.LEFT|wx.RIGHT, border=5)
		sizer.Add(skC,     pos=(row+7, 0), span=(1, 1), flag=wx.EXPAND|wx.LEFT|wx.RIGHT, border=5)
		sizer.Add(skCText, pos=(row+7, 1), span=(1, 1), flag=wx.EXPAND|wx.LEFT|wx.RIGHT, border=5)
		sizer.Add(skCDec,  pos=(row+7, 2), span=(1, 1), flag=wx.EXPAND|wx.LEFT|wx.RIGHT, border=5)
		sizer.Add(skCInc,  pos=(row+7, 3), span=(1, 1), flag=wx.EXPAND|wx.LEFT|wx.RIGHT, border=5)
		
		line = wx.StaticLine(panel)
		sizer.Add(line, pos=(row+8, 0), span=(1, 4), flag=wx.EXPAND|wx.BOTTOM, border=10)
		
		row += 9
		
		#
		# Buttons
		#
		
		ok = wx.Button(panel, ID_MASKING_OK, 'Ok', size=(56, 28))
		sizer.Add(ok, pos=(row+0, 3), flag=wx.RIGHT|wx.BOTTOM, border=5)
		
		panel.SetSizerAndFit(sizer)

		self.bpRText = bpRText
		self.dcPText = dcPText
		self.dcCText = dcCText
		self.skSText = skSText
		self.skCText = skCText
		
		#
		# Set current values
		#
		self.bpRText.SetValue('%.2f' % self.parent.data.bandpassCut)
		self.dcPText.SetValue('%i'   % self.parent.data.driftOrder)
		self.dcCText.SetValue('%i'   % self.parent.data.driftCut)
		self.skSText.SetValue('%i'   % self.parent.data.kurtosisSec)
		self.skCText.SetValue('%i'   % self.parent.data.kurtosisCut)
		
	def initEvents(self):
		self.Bind(wx.EVT_BUTTON, self.onBPRDecrease, id=ID_BANDPASS_CUT_DEC)
		self.Bind(wx.EVT_BUTTON, self.onBPRIncrease, id=ID_BANDPASS_CUT_INC)
		
		self.Bind(wx.EVT_BUTTON, self.onDCPDecrease, id=ID_DRIFT_POLY_DEC)
		self.Bind(wx.EVT_BUTTON, self.onDCPIncrease, id=ID_DRIFT_POLY_INC)
		self.Bind(wx.EVT_BUTTON, self.onDCCDecrease, id=ID_DRIFT_CUT_DEC)
		self.Bind(wx.EVT_BUTTON, self.onDCCIncrease, id=ID_DRIFT_CUT_INC)
		
		self.Bind(wx.EVT_BUTTON, self.onSKSDecrease, id=ID_SK_SEC_DEC)
		self.Bind(wx.EVT_BUTTON, self.onSKSIncrease, id=ID_SK_SEC_INC)
		self.Bind(wx.EVT_BUTTON, self.onSKCDecrease, id=ID_SK_CUT_DEC)
		self.Bind(wx.EVT_BUTTON, self.onSKCIncrease, id=ID_SK_CUT_INC)
		
		self.Bind(wx.EVT_BUTTON, self.onOk, id=ID_MASKING_OK)
		
	def onBPRDecrease(self, event):
		if self.parent.data.bandpassCut > 0.05:
			self.parent.data.bandpassCut -= 0.05
			self.bpRText.SetValue('%.2f' % self.parent.data.bandpassCut)
		
	def onBPRIncrease(self, event):
		if self.parent.data.bandpassCut < 1:
			self.parent.data.bandpassCut += 0.05
			self.bpRText.SetValue('%.2f' % self.parent.data.bandpassCut)
		
	def onDCPDecrease(self, event):
		if self.parent.data.driftOrder > 1:
			self.parent.data.driftOrder -= 1
			self.dcPText.SetValue('%i'   % self.parent.data.driftOrder)
		
	def onDCPIncrease(self, event):
		if self.parent.data.driftOrder < 12:
			self.parent.data.driftOrder += 1
			self.dcPText.SetValue('%i'   % self.parent.data.driftOrder)
	
	def onDCCDecrease(self, event):
		if self.parent.data.driftCut > 2:
			self.parent.data.driftCut -= 1
			self.dcCText.SetValue('%i'   % self.parent.data.driftCut)
		
	def onDCCIncrease(self, event):
		if self.parent.data.driftOrder < numpy.ceil(self.parent.data.data.spec.shape[0]/300):
			self.parent.data.driftCut += 1
			self.dcCText.SetValue('%i'   % self.parent.data.driftCut)
			
	def onSKSDecrease(self, event):
		if self.parent.data.kurtosisSec > 1:
			self.parent.data.kurtosisSec -= 1
			self.skSText.SetValue('%i'   % self.parent.data.kurtosisSec)
		
	def onSKSIncrease(self, event):
		if self.parent.data.kurtosisSec < 30:
			self.parent.data.kurtosisSec += 1
			self.skSText.SetValue('%i'   % self.parent.data.kurtosisSec)
	
	def onSKCDecrease(self, event):
		if self.parent.data.kurtosisCut > 2:
			self.parent.data.kurtosisCut -= 1
			self.skCText.SetValue('%i'   % self.parent.data.kurtosisCut)
		
	def onSKCIncrease(self, event):
		if self.parent.data.kurtosisCut < 12:
			self.parent.data.kurtosisCut += 1
			self.skCText.SetValue('%i'   % self.parent.data.kurtosisCut)
	
	def onOk(self, event):
		self.parent.cAdjust = None
		self.Close()


class WaterfallDisplay(wx.Frame):
	"""
	Window for displaying the waterfall data in a zoomable fashion
	"""
	
	def __init__(self, parent):
		wx.Frame.__init__(self, parent, title='Waterfall', size=(400, 375))
		
		self.parent = parent
		
		self.initUI()
		self.initEvents()
		self.Show()
		
		self.initPlot()
		
	def initUI(self):
		"""
		Start the user interface.
		"""
		
		self.statusbar = self.CreateStatusBar()
		
		hbox = wx.BoxSizer(wx.HORIZONTAL)
		
		# Add plots to panel 1
		panel1 = wx.Panel(self, -1)
		vbox1 = wx.BoxSizer(wx.VERTICAL)
		self.figure = Figure()
		self.canvas = FigureCanvasWxAgg(panel1, -1, self.figure)
		self.toolbar = NavigationToolbar2WxAgg(self.canvas)
		self.toolbar.Realize()
		vbox1.Add(self.canvas,  1, wx.EXPAND)
		vbox1.Add(self.toolbar, 0, wx.LEFT | wx.FIXED_MINSIZE)
		panel1.SetSizer(vbox1)
		hbox.Add(panel1, 1, wx.EXPAND)
		
		# Use some sizers to see layout options
		self.SetSizer(hbox)
		self.SetAutoLayout(1)
		hbox.Fit(self)
		
	def initEvents(self):
		"""
		Set all of the various events in the data range window.
		"""
		
		# Make the images resizable
		self.Bind(wx.EVT_PAINT, self.resizePlots)
		
	def initPlot(self):
		"""
		Populate the figure/canvas areas with a plot.  We only need to do this
		once for this type of window.
		"""

		if self.parent.data.index / (self.parent.data.spec.shape[1]/2) == 0:
			freq = self.parent.data.freq1
		else:
			freq = self.parent.data.freq2
		
		if self.parent.data.bandpass:
			spec = self.parent.data.specBandpass
			limits = self.parent.data.limitsBandpass
		else:
			spec = self.parent.data.spec
			limits = self.parent.data.limits
		
		# Plot Waterfall
		self.figure.clf()
		self.ax1 = self.figure.gca()
		if self.parent.data.usedB:
			m = self.ax1.imshow(to_dB(spec[:,self.parent.data.index,:]), interpolation='nearest', extent=(freq[0]/1e6, freq[-1]/1e6, self.parent.data.time[0], self.parent.data.time[-1]), origin='lower', vmin=limits[self.parent.data.index][0], vmax=limits[self.parent.data.index][1])
			try:
				cm = self.figure.colorbar(m, use_gridspec=True)
			except:
				cm = self.figure.colorbar(m)
			cm.ax.set_ylabel('PSD [arb. dB]')
		else:
			m = self.ax1.imshow(spec[:,self.parent.data.index,:], interpolation='nearest', extent=(freq[0]/1e6, freq[-1]/1e6, self.parent.data.time[0], self.parent.data.time[-1]), origin='lower', vmin=limits[self.parent.data.index][0], vmax=limits[self.parent.data.index][1])
			try:
				cm = self.figure.colorbar(m, use_gridspec=True)
			except:
				cm = self.figure.colorbar(m)
			cm.ax.set_ylabel('PSD [arb. lin.]')
		self.ax1.axis('auto')
		self.ax1.set_xlim((freq[0]/1e6, freq[-1]/1e6))
		self.ax1.set_ylim((self.parent.data.time[0], self.parent.data.time[-1]))
		self.ax1.set_xlabel('Frequency [MHz]')
		self.ax1.set_ylabel('Elapsed Time - %.3f [s]' % (self.parent.data.iOffset*self.parent.data.tInt))
		if self.parent.data.linear:
			tun = self.parent.data.index / 4 + 1
			ind = self.parent.data.index % 4
			mapper = {0: 'XX', 1: 'XY', 2: 'YX', 3: 'YY'}
			self.ax1.set_title('Tuning %i, %s' % (tun, mapper[ind]))
		else:
			tun = self.parent.data.index / 4 + 1
			ind = self.parent.data.index % 4
			mapper = {0: 'I', 1: 'Q', 2: 'U', 3: 'V'}
			self.ax1.set_title('Tuning %i, %s' % (tun, mapper[ind]))
			
		## Draw and save the click (Why?)
		self.canvas.draw()
		self.connect()
		
	def connect(self):
		"""
		Connect to all the events we need to interact with the plots.
		"""
		
		self.cidmotion  = self.figure.canvas.mpl_connect('motion_notify_event', self.on_motion)
	
	def on_motion(self, event):
		"""
		Deal with motion events in the stand field window.  This involves 
		setting the status bar with the current x and y coordinates as well
		as the stand number of the selected stand (if any).
		"""

		if self.parent.data.index / (self.parent.data.spec.shape[1]/2) == 0:
			freq = self.parent.data.freq1
		else:
			freq = self.parent.data.freq2
		
		if event.inaxes:
			clickX = event.xdata
			clickY = event.ydata
			
			dataX = numpy.where(numpy.abs(clickX-freq/1e6) == (numpy.abs(clickX-freq/1e6).min()))[0][0]
			dataY = numpy.where(numpy.abs(clickY-self.parent.data.time) == (numpy.abs(clickY-self.parent.data.time).min()))[0][0]
			
			if self.parent.data.usedB:
				value = to_dB(self.parent.data.spec[dataY, self.parent.data.index, dataX])
				self.statusbar.SetStatusText("f=%.4f MHz, t=%.4f s, p=%.2f dB" % (clickX, clickY, value))
			else:
				value = self.parent.data.spec[dataY, self.parent.data.index, dataX]
				self.statusbar.SetStatusText("f=%.4f MHz, t=%.4f s, p=%.2f" % (clickX, clickY, value))
		else:
			self.statusbar.SetStatusText("")
	
	def disconnect(self):
		"""
		Disconnect all the stored connection ids.
		"""
		
		self.figure.canvas.mpl_disconnect(self.cidmotion)
		
	def onCancel(self, event):
		self.Close()
		
	def resizePlots(self, event):
		w, h = self.GetSize()
		dpi = self.figure.get_dpi()
		newW = 1.0*w/dpi
		newH1 = 1.0*(h/2-100)/dpi
		newH2 = 1.0*(h/2-75)/dpi
		self.figure.set_size_inches((newW, newH1))
		self.figure.canvas.draw()

	def GetToolBar(self):
		# You will need to override GetToolBar if you are using an 
		# unmanaged toolbar in your frame
		return self.toolbar


class DriftCurveDisplay(wx.Frame):
	"""
	Window for displaying the waterfall data in a zoomable fashion
	"""
	
	def __init__(self, parent):
		wx.Frame.__init__(self, parent, title='Waterfall', size=(400, 375))
		
		self.parent = parent
		
		self.initUI()
		self.initEvents()
		self.Show()
		
		self.initPlot()
		
		self.site = stations.lwa1.getObserver()
		
	def initUI(self):
		"""
		Start the user interface.
		"""
		
		self.statusbar = self.CreateStatusBar()
		
		hbox = wx.BoxSizer(wx.HORIZONTAL)
		
		# Add plots to panel 1
		panel1 = wx.Panel(self, -1)
		vbox1 = wx.BoxSizer(wx.VERTICAL)
		self.figure = Figure()
		self.canvas = FigureCanvasWxAgg(panel1, -1, self.figure)
		self.toolbar = NavigationToolbar2WxAgg(self.canvas)
		self.toolbar.Realize()
		vbox1.Add(self.canvas,  1, wx.EXPAND)
		vbox1.Add(self.toolbar, 0, wx.LEFT | wx.FIXED_MINSIZE)
		panel1.SetSizer(vbox1)
		hbox.Add(panel1, 1, wx.EXPAND)
		
		# Use some sizers to see layout options
		self.SetSizer(hbox)
		self.SetAutoLayout(1)
		hbox.Fit(self)
		
	def initEvents(self):
		"""
		Set all of the various events in the data range window.
		"""
		
		# Make the images resizable
		self.Bind(wx.EVT_PAINT, self.resizePlots)
		
	def initPlot(self):
		"""
		Populate the figure/canvas areas with a plot.  We only need to do this
		once for this type of window.
		"""
		
		if self.parent.data.bandpass:
			spec = self.parent.data.specBandpass
			limits = self.parent.data.limitsBandpass
		else:
			spec = self.parent.data.spec
			limits = self.parent.data.limits
		
		# Plot Drift Curve
		self.figure.clf()
		self.ax1 = self.figure.gca()
		
		self.drift = spec[:,:,spec.shape[2]/8:7*spec.shape[2]/8].sum(axis=2)
		
		if self.parent.data.usedB:
			self.ax1.plot(self.parent.data.time, to_dB(self.drift[:,self.parent.data.index]), linestyle='-', marker='x')
			self.ax1.set_ylabel('Inner 75% Total Power [arb. dB]')
		else:
			self.ax1.plot(self.parent.data.time, self.drift[:,self.parent.data.index], linestyle='-', marker='x')
			self.ax1.set_ylabel('Inner 75% Total Power [arb. lin.]')
		self.ax1.set_xlim((self.parent.data.time[0], self.parent.data.time[-1]))
		self.ax1.set_ylabel('Elapsed Time - %.3f [s]' % (self.parent.data.iOffset*self.parent.data.tInt))
		if self.parent.data.linear:
			tun = self.parent.data.index / 4 + 1
			ind = self.parent.data.index % 4
			mapper = {0: 'XX', 1: 'XY', 2: 'YX', 3: 'YY'}
			self.ax1.set_title('Tuning %i, %s' % (tun, mapper[ind]))
		else:
			tun = self.parent.data.index / 4 + 1
			ind = self.parent.data.index % 4
			mapper = {0: 'I', 1: 'Q', 2: 'U', 3: 'V'}
			self.ax1.set_title('Tuning %i, %s' % (tun, mapper[ind]))
			
		## Draw and save the click (Why?)
		self.canvas.draw()
		self.connect()
		
	def connect(self):
		"""
		Connect to all the events we need to interact with the plots.
		"""
		
		self.cidmotion  = self.figure.canvas.mpl_connect('motion_notify_event', self.on_motion)
	
	def on_motion(self, event):
		"""
		Deal with motion events in the stand field window.  This involves 
		setting the status bar with the current x and y coordinates as well
		as the stand number of the selected stand (if any).
		"""
		
		if event.inaxes:
			clickX = event.xdata
			clickY = event.ydata
			
			dataX = numpy.where(numpy.abs(clickX-self.parent.data.time) == (numpy.abs(clickX-self.parent.data.time).min()))[0][0]
			
			ts = datetime.utcfromtimestamp(self.parent.data.timesNPZ[dataX])
			self.site.date = ts.strftime('%Y/%m/%d %H:%M:%S')
			lst = self.site.sidereal_time()
			
			if self.parent.data.usedB:
				value = to_dB(self.drift[dataX,self.parent.data.index])
				self.statusbar.SetStatusText("t=%s, LST=%s, p=%.2f dB" % (ts, lst, value))
			else:
				value = self.drift[dataX,self.parent.data.index]
				self.statusbar.SetStatusText("t=%s, LST=%s, p=%.2f" % (ts, lst, value))
		else:
			self.statusbar.SetStatusText("")
	
	def disconnect(self):
		"""
		Disconnect all the stored connection ids.
		"""
		
		self.figure.canvas.mpl_disconnect(self.cidmotion)
		
	def onCancel(self, event):
		self.Close()
		
	def resizePlots(self, event):
		w, h = self.GetSize()
		dpi = self.figure.get_dpi()
		newW = 1.0*w/dpi
		newH1 = 1.0*(h/2-100)/dpi
		newH2 = 1.0*(h/2-75)/dpi
		self.figure.set_size_inches((newW, newH1))
		self.figure.canvas.draw()

	def GetToolBar(self):
		# You will need to override GetToolBar if you are using an 
		# unmanaged toolbar in your frame
		return self.toolbar


def main(args):
	# Turn off all NumPy warnings to keep stdout clean
	errs = numpy.geterr()
	numpy.seterr(invalid='ignore', divide='ignore')
	
	# Parse the command line options
	config = parseOptions(args)
	
	# Go!
	app = wx.App(0)
	frame = MainWindow(None, -1)
	if len(config['args']) == 1:
		frame.filename = config['args'][0]
		frame.offset = config['offset']
		frame.duration = config['duration']
		frame.data = Waterfall_GUI(frame)
		frame.data.loadData(config['args'][0], obsID=config['obsID'])
		frame.render()
		frame.data.render()
		frame.data.draw()
		
		if frame.data.filenames is None:
			frame.examineFileButton.Enable(False)
		else:
			frame.examineFileButton.Enable(True)
	app.MainLoop()


if __name__ == "__main__":
	main(sys.argv[1:])
