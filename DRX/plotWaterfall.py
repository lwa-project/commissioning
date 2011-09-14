#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Given the NPZ output of drxWaterfall, plot it in an interative way.

$Rev$
$LastChangedBy$
$LastChangedDate$
"""

import os
import sys
import time
import numpy
from datetime import datetime

from lsl.common.dp import fS
from lsl.common import stations
from lsl.misc.mathutil import to_dB, from_dB, savitzky_golay
from lsl.statistics import robust

import wx
import matplotlib
matplotlib.use('WXAgg')
matplotlib.interactive(True)

from matplotlib.backends.backend_wxagg import NavigationToolbar2WxAgg, FigureCanvasWxAgg
from matplotlib.figure import Figure


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


def skStd(M, N=1):
	"""
	Return the expected standard deviation of the spectral kurtosis for M points 
	each composed of N measurments.
	
	.. note::
		In the future (LSL 0.5), this will be lsl.statistics.kurtosis.std()
	"""
	
	return numpy.sqrt( skVar(M, N) )


def skVar(M, N=1):
	"""
	Return the expected variance (second central moment) of the spectral kurtosis 
	for M points each composed of N measurments.
	.. note::
		In the future (LSL 0.5), this will be lsl.statistics.kurtosis.var()
	"""

	return 2.0*N*(N+1)*M**2/ float( (M-1)*(M*N+3)*(M*N+2) )


class Waterfall_GUI(object):
	def __init__(self, frame, freq=None, spec=None, tInt=None):
		self.frame = frame
		self.press = None
		
		self.filename = ''
		self.index = 0
		
		self.bandpass = False
		self.freq = freq
		self.spec = spec
		self.tInt = None
		
		self.ax1a = None
		self.ax1b = None
		self.ax2 = None
		
		self.spectrumClick = None
		
	def loadData(self, filename):
		"""
		Load in data from an NPZ file.
		"""
		
		print "Loading file '%s'" % os.path.split(filename)[1]
		tStart = time.time()
		
		# Save the filename
		self.filename = filename
		dataDictionary = numpy.load(filename)
		
		# Load the Data
		print " %6.3f s - Extracting data" % (time.time() - tStart)
		try:
			self.srate = dataDictionary['srate']
		except KeyError:
			self.srate = 19.6e6
		self.freq = dataDictionary['freq']
		try:
			mask = dataDictionary['mask'].astype(numpy.bool)
		except KeyError:
			mask = numpy.zeros(dataDictionary['spec'].shape, dtype=numpy.bool)
		self.spec = numpy.ma.array(dataDictionary['spec'], mask=mask)
		self.tInt = dataDictionary['tInt']
		self.time = dataDictionary['times']
		self.time -= self.time.min()
		
		# Construct frequency and time master masks to prevent some masked things from getting unmasked
		self.freqMask = numpy.median(self.spec.mask, axis=0)
		self.timeMask = numpy.median(self.spec.mask, axis=2)
		
		# Other data to keep around in case we save
		self.timesNPZ = dataDictionary['times']
		self.standMapperNPZ = dataDictionary['standMapper']
		
		# Deal with the potential for aggregated files
		try:
			self.tIntActual = dataDictionary['tIntActual']
			self.tIntOriginal = dataDictionary['tIntOriginal']
			self.filenames = dataDictionary['filenames']
			print "*** NPZ file appears to be aggregation of %i files ***" % len(self.filenames)
		except KeyError:
			self.tIntActual = self.tInt
			self.tIntOriginal = self.tInt
			self.filenames = None
		
		## Get the filter model
		#print " %6.3f s - Building DRX bandpass model" % (time.time() - tStart)
		#self.bpm = drxFilter(sampleRate=self.srate)(self.freq)
		
		# Compute the bandpass fit
		print " %6.3f s - Computing bandpass fits" % (time.time() - tStart)
		self.computeBandpass()
		
		# Find the mean spectra
		print " %6.3f s - Computing mean spectra" % (time.time() - tStart)
		self.mean = numpy.mean(self.spec, axis=0)
		self.meanBandpass = numpy.mean(self.specBandpass, axis=0)
		
		# Set default colobars
		print " %6.3f s - Setting default colorbar ranges" % (time.time() - tStart)
		self.limits = []
		for i in xrange(self.spec.shape[1]):
			self.limits.append( [to_dB(self.spec[:,i,:]).min(), to_dB(self.spec[:,i,:]).max()] )
			
		self.limitsBandpass = []
		toUse = numpy.arange(self.spec.shape[2]/10, 9*self.spec.shape[2]/10)
		for i in xrange(self.spec.shape[1]):
			self.limitsBandpass.append( [to_dB(self.specBandpass[:,i,toUse]).min(), to_dB(self.specBandpass[:,i,toUse]).max()] )
		
		try:
			self.disconnect()
		except:
			pass
		
		# Clear the old marks
		self.oldMarkA = None
		self.oldMarkB = None
		self.frame.figure1a.clf()
		self.frame.figure1a.clf()
		self.frame.figure2.clf()
		
		self.connect()
		
		print " %6.3f s - Ready" % (time.time() - tStart)
	
	def computeBandpass(self):
		"""
		Compute the bandpass fits.
		"""
		
		freq2 = self.freq / 1e6
		meanSpec = numpy.mean(self.spec.data, axis=0)
		bpm2 = []
		for i in xrange(self.spec.shape[1]):
			from scipy.interpolate import splrep, splev
			#t = splrep(freq2, numpy.log10(meanSpec[i,:])*10, s=1e-5, k=3)
			#noiseSlope = 10.0**(splev(freq2, t, der=0)/10.0)
			noiseSlope = savitzky_golay(to_dB(meanSpec[i,:]), 41, 9, deriv=0)
			noiseSlope = from_dB(noiseSlope)
			
			bpm2.append( noiseSlope/noiseSlope.mean() )
			
		# Apply the bandpass correction
		self.specBandpass = numpy.ma.array(self.spec.data*1.0, mask=self.spec.mask)
		for i in xrange(self.spec.shape[1]):
			for j in xrange(self.spec.shape[0]):
				self.specBandpass[j,i,:] = self.spec[j,i,:] / bpm2[i]
	
	def draw(self):
		"""
		Draw the waterfall diagram and the total power with time.
		"""
		
		if self.bandpass:
			spec = self.specBandpass
			limits = self.limitsBandpass
		else:
			spec = self.spec
			limits = self.limits
		
		# Plot 1(a) - Waterfall
		self.frame.figure1a.clf()
		self.ax1a = self.frame.figure1a.gca()
		m = self.ax1a.imshow(to_dB(spec[:,self.index,:]), interpolation='nearest', extent=(self.freq[0]/1e6, self.freq[-1]/1e6, self.time[0], self.time[-1]), origin='lower', vmin=limits[self.index][0], vmax=limits[self.index][1])
		cm = self.frame.figure1a.colorbar(m, ax=self.ax1a)
		cm.ax.set_ylabel('PSD [arb. dB]')
		self.ax1a.axis('auto')
		self.ax1a.set_xlabel('Frequency [MHz]')
		self.ax1a.set_ylabel('Elapsed Time [s]')
		self.ax1a.set_title('Tuning %i, Pol. %s' % (self.index/2+1, 'Y' if self.index %2 else 'X'))
		
		if self.oldMarkA is not None:
			self.ax1a.lines.extend(self.oldMarkA)
		
		self.frame.canvas1a.draw()
		
		# Plot 1(b) - Drift
		self.drift = spec[:,:,spec.shape[2]/4:3*spec.shape[2]/4].sum(axis=2)
		
		self.frame.figure1b.clf()
		self.ax1b = self.frame.figure1b.gca()
		self.ax1b.plot(self.drift[:,self.index], self.time, linestyle=' ', marker='x')
		self.ax1b.set_ylim([self.time[0], self.time[-1]])
		self.ax1b.set_xlabel('Total Power [arb. dB]')
		self.ax1b.set_ylabel('Elapsed Time [s]')
		
		if self.oldMarkB is not None:
			self.ax1b.lines.extend(self.oldMarkB)
		
		self.frame.canvas1b.draw()
	
	def drawSpectrum(self, clickY, Home=False):
		"""Get the spectrum at a particular point in time."""
		
		try:
			dataY = int(round(clickY / self.tInt))
		except TypeError:
			return False
		
		if self.bandpass:
			spec = self.specBandpass[dataY,self.index,:]
			medianSpec = self.meanBandpass[self.index,:]
			limits = self.limitsBandpass
		else:
			spec = self.spec[dataY,self.index,:]
			medianSpec = self.mean[self.index,:]
			limits = self.limits
		
		if self.frame.toolbar.mode == 'zoom rect' and not Home:
			try:
				oldXlim = self.ax2.get_xlim()
				oldYlim = self.ax2.get_ylim()
			except:
				oldXlim = [self.freq[0]/1e6, self.freq[-1]/1e6]
				oldYlim = limits[self.index]
		else:
			oldXlim = [self.freq[0]/1e6, self.freq[-1]/1e6]
			oldYlim = limits[self.index]
		
		self.frame.figure2.clf()
		self.ax2 = self.frame.figure2.gca()
		self.ax2.plot(self.freq/1e6, to_dB(spec), linestyle=' ', marker='o', label='Current', mec='blue', mfc='None')
		self.ax2.plot(self.freq/1e6, to_dB(medianSpec), label='Mean', alpha=0.5, color='green')
		self.ax2.set_xlim(oldXlim)
		self.ax2.set_ylim(oldYlim)
		self.ax2.legend(loc=0)
		self.ax2.set_xlabel('Frequency [MHz]')
		self.ax2.set_ylabel('PSD [arb. dB]')
		
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
		
		oldXSizeA = self.ax1a.get_xlim()
		oldYSizeA = self.ax1a.get_ylim()
		
		oldXSizeB = self.ax1b.get_xlim()
		oldYSizeB = self.ax1b.get_ylim()
		
		self.oldMarkA = self.ax1a.plot([-1e6, 1e6], [self.time[dataY],]*2, color='red')
		self.oldMarkB = self.ax1b.plot([-1e6, 1e6], [self.time[dataY],]*2, color='red')
		
		self.ax1a.set_xlim(oldXSizeA)
		self.ax1a.set_ylim(oldYSizeA)
		
		self.ax1b.set_xlim(oldXSizeB)
		self.ax1b.set_ylim(oldYSizeB)
		
		self.frame.canvas1a.draw()
		self.frame.canvas1b.draw()
	
	def connect(self):
		"""
		Connect to all the events we need
		"""
		
		self.cidpress1a = self.frame.figure1a.canvas.mpl_connect('button_press_event', self.on_press1a)
		self.cidpress1b = self.frame.figure1b.canvas.mpl_connect('button_press_event', self.on_press1b)
		self.cidpress2  = self.frame.figure2.canvas.mpl_connect('button_press_event', self.on_press2)
		self.cidmotion  = self.frame.figure1a.canvas.mpl_connect('motion_notify_event', self.on_motion)
		self.frame.toolbar.home = self.on_home
		
	def on_home(self, *args):
		"""
		Override of the toolbar 'home' operation.
		"""
		
		self.drawSpectrum(self.spectrumClick, Home=True)
	
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
			
			scaleX = self.ax1b.get_xlim()
			rangeX = scaleX[1] - scaleX[0]
			
			scaleY = self.ax1b.get_ylim()
			rangeY = scaleY[1] - scaleY[0]
			
			dataY = int(round(clickY / self.tInt))
			lower = dataY - 200
			lower = 0 if lower < 0 else lower
			upper = dataY + 200
			upper = self.drift.shape[0]-1 if upper > self.drift.shape[0]-1 else upper
			
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
		if event.inaxes:
			clickX = event.xdata
			clickY = event.ydata
			
			dataX = numpy.where(numpy.abs(clickX-self.freq/1e6) == (numpy.abs(clickX-self.freq/1e6).min()))[0][0]
			
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
		if event.inaxes:
			clickX = event.xdata
			clickY = event.ydata
			
			dataX = numpy.where(numpy.abs(clickX-self.freq/1e6) == (numpy.abs(clickX-self.freq/1e6).min()))[0][0]
			dataY = numpy.where(numpy.abs(clickY-self.time) == (numpy.abs(clickY-self.time).min()))[0][0]
			
			value = to_dB(self.spec[dataY, self.index, dataX])
			self.frame.statusbar.SetStatusText("f=%.4f MHz, t=%.4f s, p=%.2f dB" % (clickX, clickY, value))
		else:
			self.frame.statusbar.SetStatusText("")
			
	
	def disconnect(self):
		'disconnect all the stored connection ids'
		
		self.frame.figure1a.canvas.mpl_disconnect(self.cidpress1a)
		self.frame.figure1b.canvas.mpl_disconnect(self.cidpress1b)
		self.frame.figure2.canvas.mpl_disconnect(self.cidpress2)
		self.frame.figure1a.canvas.mpl_disconnect(self.cidmotion)


ID_OPEN    = 10
ID_SAVE    = 11
ID_SAVE_AS = 12
ID_QUIT    = 13
ID_COLOR_AUTO = 20
ID_COLOR_ADJUST = 21
ID_TUNING1_X = 30
ID_TUNING1_Y = 31
ID_TUNING2_X = 32
ID_TUNING2_Y = 33
ID_ZOOMABLE_WATERFALL = 34
ID_ZOOMABLE_DRIFT = 35
ID_MASK_SUGGEST_CURRENT = 40
ID_MASK_SUGGEST_ALL = 41
ID_MASK_RESET_CURRENT = 42
ID_MASK_RESET_ALL = 43
ID_MASK_TWEAK = 44
ID_BANDPASS_ON = 50
ID_BANDPASS_OFF = 51
ID_BANDPASS_RECOMPUTE = 52

class MainWindow(wx.Frame):
	def __init__(self, parent, id):
		self.dirname = ''
		self.filename = ''
		self.data = None
		
		wx.Frame.__init__(self, parent, id, title="DRX Waterfall Viewer", size=(600,800))
		
		self.initUI()
		self.initEvents()
		self.Show()
		
		self.bandpassCut = 0.85
		self.driftOrder  = 5
		self.driftCut    = 4
		self.kurtosisSec = 1
		self.kurtosisCut = 4
		
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
		
		## File Menu
		open = wx.MenuItem(fileMenu, ID_OPEN, "&Open")
		fileMenu.AppendItem(open)
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
		dataMenu.AppendRadioItem(ID_TUNING1_X, 'Tuning 1, Pol. X')
		dataMenu.AppendRadioItem(ID_TUNING1_Y, 'Tuning 1, Pol. Y')
		dataMenu.AppendSeparator()
		dataMenu.AppendRadioItem(ID_TUNING2_X, 'Tuning 2, Pol. X')
		dataMenu.AppendRadioItem(ID_TUNING2_Y, 'Tuning 2, Pol. Y')
		dataMenu.AppendSeparator()
		zm = wx.MenuItem(dataMenu, ID_ZOOMABLE_WATERFALL, 'Zoomable Waterfall')
		dataMenu.AppendItem(zm)
		zd = wx.MenuItem(dataMenu, ID_ZOOMABLE_DRIFT, 'Zoomable Drift Curve')
		dataMenu.AppendItem(zd)
		
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

		# Creating the menubar.
		menuBar.Append(fileMenu,"&File") # Adding the "filemenu" to the MenuBar
		menuBar.Append(colorMenu, "&Color")
		menuBar.Append(dataMenu, "&Data")
		menuBar.Append(maskMenu, "&Mask")
		menuBar.Append(bandpassMenu, "&Bandpass")
		self.SetMenuBar(menuBar)  # Adding the MenuBar to the Frame content.
		
		vbox = wx.BoxSizer(wx.VERTICAL)
		
		# Add waterfall plot
		panel1 = wx.Panel(self, -1)
		hbox1 = wx.BoxSizer(wx.HORIZONTAL)
		self.figure1a = Figure()
		self.canvas1a = FigureCanvasWxAgg(panel1, -1, self.figure1a)
		hbox1.Add(self.canvas1a, 1, wx.EXPAND)
		
		# Add a total power with time plot
		self.figure1b = Figure()
		self.canvas1b = FigureCanvasWxAgg(panel1, -1, self.figure1b)
		hbox1.Add(self.canvas1b, 1, wx.EXPAND)
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
		hbox3.Add(self.toolbar, 0, wx.LEFT | wx.FIXED_MINSIZE)
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
		
		self.Bind(wx.EVT_MENU, self.onTuning1X, id=ID_TUNING1_X)
		self.Bind(wx.EVT_MENU, self.onTuning1Y, id=ID_TUNING1_Y)
		self.Bind(wx.EVT_MENU, self.onTuning2X, id=ID_TUNING2_X)
		self.Bind(wx.EVT_MENU, self.onTuning2Y, id=ID_TUNING2_Y)
		self.Bind(wx.EVT_MENU, self.onZoomWaterfall, id=ID_ZOOMABLE_WATERFALL)
		self.Bind(wx.EVT_MENU, self.onZoomDrift, id=ID_ZOOMABLE_DRIFT)
		
		self.Bind(wx.EVT_MENU, self.onMaskSuggestCurrent, id=ID_MASK_SUGGEST_CURRENT)
		self.Bind(wx.EVT_MENU, self.onMaskSuggestAll, id=ID_MASK_SUGGEST_ALL)
		self.Bind(wx.EVT_MENU, self.onMaskResetCurrent, id=ID_MASK_RESET_CURRENT)
		self.Bind(wx.EVT_MENU, self.onMaskResetAll, id=ID_MASK_RESET_ALL)
		self.Bind(wx.EVT_MENU, self.onMaskTweak, id=ID_MASK_TWEAK)
		
		self.Bind(wx.EVT_MENU, self.onBandpassOn, id=ID_BANDPASS_ON)
		self.Bind(wx.EVT_MENU, self.onBandpassOff, id=ID_BANDPASS_OFF)
		self.Bind(wx.EVT_MENU, self.onBandpassRecompute, id=ID_BANDPASS_RECOMPUTE)
		
		# Key events
		self.canvas1a.Bind(wx.EVT_KEY_UP, self.onKeyPress)
		self.canvas1b.Bind(wx.EVT_KEY_UP, self.onKeyPress)
		self.canvas2.Bind(wx.EVT_KEY_UP,  self.onKeyPress)
		
		# Make the images resizable
		self.Bind(wx.EVT_PAINT, self.resizePlots)
	
	def onOpen(self, event):
		"""
		Open a file.
		"""
		
		dlg = wx.FileDialog(self, "Choose a file", self.dirname, "", "*.*", wx.OPEN)
		if dlg.ShowModal() == wx.ID_OK:
			self.filename = dlg.GetFilename()
			self.dirname = dlg.GetDirectory()
			self.data = Waterfall_GUI(self)
			self.data.loadData(os.path.join(self.dirname, self.filename))
			self.data.draw()
			
			if self.cAdjust is not None:
				try:
					self.cAdjust.Close()
				except:
					pass
				self.cAdjust = None
		dlg.Destroy()
		
	def onSave(self, event):
		"""
		Save the data mask to a new NPZ file.
		"""
		
		if self.data.filename == '':
			self.onSaveAs(event)
		else:
			numpy.savez(self.filename, freq=self.data.freq, times=self.data.timesNPZ, spec=self.data.spec.data, mask=self.data.spec.mask, tInt=self.data.tInt,  standMapper=self.data.standMapperNPZ)

	def onSaveAs(self, event):
		"""
		Save the current observation to a new SD file.
		"""
		
		dialog = wx.FileDialog(self, "Select Output File", self.dirname, '', 'NPZ Files (*.npz)|*.npz|All Files (*.*)|*.*', wx.SAVE|wx.FD_OVERWRITE_PROMPT)
			
		if dialog.ShowModal() == wx.ID_OK:
			self.dirname = dialog.GetDirectory()
			self.filename = dialog.GetPath()
			
			numpy.savez(os.path.join(self.dirname, self.filename), freq=self.data.freq, times=self.data.timesNPZ, spec=self.data.spec.data, mask=self.data.spec.mask, tInt=self.data.tInt,  standMapper=self.data.standMapperNPZ)
			
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
		
		i = self.data.index
		toUse = numpy.arange(self.data.spec.shape[2]/10, 9*self.data.spec.shape[2]/10)
		if self.data.bandpass:
			self.data.limitsBandpass[i] = [to_dB(self.data.specBandpass[:,i,toUse]).min(), to_dB(self.data.specBandpass[:,i,toUse]).max()] 
		else:
			self.data.limits[i] = [to_dB(self.data.spec[:,i,:]).min(), to_dB(self.data.spec[:,i,:]).max()]
			
		self.data.draw()
		self.data.drawSpectrum(self.data.spectrumClick)
		self.data.makeMark(self.data.spectrumClick)
		
		wx.EndBusyCursor()
		
	def onAdjust(self, event):
		"""
		Bring up the colorbar adjustment dialog window.
		"""
		
		ContrastAdjust(self)
		
	def onTuning1X(self, event):
		"""
		Display tuning 1, pol X.
		"""
		
		wx.BeginBusyCursor()
		
		self.data.index = 0
		self.data.draw()
		if self.data.spectrumClick is not None:
			self.data.drawSpectrum(self.data.spectrumClick)
			
		wx.EndBusyCursor()
		
	def onTuning1Y(self, event):
		"""
		Display tuning 1, pol Y.
		"""
		
		wx.BeginBusyCursor()
		
		self.data.index = 1
		self.data.draw()
		if self.data.spectrumClick is not None:
			self.data.drawSpectrum(self.data.spectrumClick)
			
		wx.EndBusyCursor()
		
	def onTuning2X(self, event):
		"""
		Display tuning 2, pol X.
		"""
		
		wx.BeginBusyCursor()
		
		self.data.index = 2
		self.data.draw()
		if self.data.spectrumClick is not None:
			self.data.drawSpectrum(self.data.spectrumClick)
			
		wx.EndBusyCursor()
		
	def onTuning2Y(self, event):
		"""
		Display tuning 2, pol Y.
		"""
		
		wx.BeginBusyCursor()
		
		self.data.index = 3
		self.data.draw()
		if self.data.spectrumClick is not None:
			self.data.drawSpectrum(self.data.spectrumClick)
			
		wx.EndBusyCursor()
		
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
		
	def onMaskSuggestCurrent(self, event):
		"""
		Suggest a series of frequency and time-based masks to apply to 
		the current tuning/polarization.
		"""
		
		wx.BeginBusyCursor()
		
		bad = numpy.where( (self.data.freq < -self.data.srate/2*self.bandpassCut) | (self.data.freq > self.data.srate/2*self.bandpassCut) )[0]
		for b in bad:
			self.data.spec.mask[:,self.data.index,b] = True
			self.data.specBandpass.mask[:,self.data.index,b] = True
			self.data.freqMask[self.data.index,b] = True
		
		spec = self.data.spec.data[:,self.data.index,:]
		drift = spec.sum(axis=1)
		coeff = robust.polyfit(numpy.arange(drift.size), drift, self.driftOrder)
		fit = numpy.polyval(coeff, numpy.arange(drift.size))
		rDrift = drift / fit
		
		mean = robust.mean(rDrift)
		std  = robust.std(rDrift)
		
		bad = numpy.where( numpy.abs(rDrift - mean) >= self.driftCut*std )[0]
		for b in bad:
			self.data.spec.mask[b,self.data.index,:] = True
			self.data.specBandpass.mask[b,self.data.index,:] = True
			self.data.timeMask[b,self.data.index] = True
		
		N = self.data.srate/(self.data.freq.size+1)*self.data.tIntActual
		kurtosis = numpy.zeros((self.kurtosisSec, self.data.spec.shape[2]))
		
		secSize = self.data.spec.shape[0]/self.kurtosisSec
		for k in xrange(self.kurtosisSec):
			tStart = k*secSize
			tStop  = (k+1)*secSize
			
			for j in xrange(self.data.spec.shape[2]):
				channel = self.data.spec.data[tStart:tStop,self.data.index,j]
				kurtosis[k,j] = spectralKurtosis(channel, N=N)
		
		kMean = 1.0
		kStd  = skStd(secSize, N)
		print self.data.tIntActual, kStd, kurtosis.min(), kurtosis.max()
		
		bad = numpy.where( numpy.abs(kurtosis - kMean) >= self.kurtosisCut*kStd )
		for k,b in zip(bad[0], bad[1]):
			tStart = k*secSize
			tStop  = (k+1)*secSize
			
			try:
				for j in xrange(b-2, b+3):
					self.data.spec.mask[tStart:tStop,self.data.index,j] = True
					self.data.specBandpass.mask[tStart:tStop,self.data.index,j] = True
					self.data.freqMask[self.data.index,j] = True
			except IndexError:
				pass
			
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
		
		for i in xrange(self.data.spec.shape[1]):
			bad = numpy.where( (self.data.freq < -self.data.srate/2*self.bandpassCut) | (self.data.freq > self.data.srate/2*self.bandpassCut) )[0]
			for b in bad:
				self.data.spec.mask[:,i,b] = True
				self.data.specBandpass.mask[:,i,b] = True
				self.data.freqMask[i,b] = True
			
			spec = self.data.spec.data[:,i,:]
			drift = spec.sum(axis=1)
			coeff = robust.polyfit(numpy.arange(drift.size), drift, self.driftOrder)
			fit = numpy.polyval(coeff, numpy.arange(drift.size))
			rDrift = drift / fit
			
			mean = robust.mean(rDrift)
			std  = robust.std(rDrift)
			
			bad = numpy.where( numpy.abs(rDrift - mean) >= self.driftCut*std )[0]
			for b in bad:
				self.data.spec.mask[b,i,:] = True
				self.data.specBandpass.mask[b,i,:] = True
				self.data.timeMask[b,i] = True
			
			N = self.data.srate/(self.data.freq.size+1)*self.data.tIntActual
			kurtosis = numpy.zeros((self.kurtosisSec, self.data.spec.shape[2]))
			
			secSize = self.data.spec.shape[0]/self.kurtosisSec
			for k in xrange(self.kurtosisSec):
				tStart = k*secSize
				tStop  = (k+1)*secSize
				
				for j in xrange(self.data.spec.shape[2]):
					channel = self.data.spec.data[tStart:tStop,i,j]
					kurtosis[k,j] = spectralKurtosis(channel, N=N)
			
			kMean = 1.0
			kStd  = skStd(secSize, N)
			
			bad = numpy.where( numpy.abs(kurtosis - kMean) >= self.kurtosisCut*kStd )
			for k,b in zip(bad[0], bad[1]):
				tStart = k*secSize
				tStop  = (k+1)*secSize
				
				try:
					for j in xrange(b-2, b+3):
						self.data.spec.mask[tStart:tStop,i,j] = True
						self.data.specBandpass.mask[tStart:tStop,i,j] = True
						self.data.freqMask[i,j] = True
				except IndexError:
					pass
			
		self.data.draw()
		self.data.drawSpectrum(self.data.spectrumClick)
		self.data.makeMark(self.data.spectrumClick)
		
		wx.EndBusyCursor()
	
	def onMaskResetCurrent(self, event):
		"""
		Reset the current mask.
		"""
		
		wx.BeginBusyCursor()
		
		self.data.spec.mask[:,self.data.index,:] = False
		self.data.specBandpass.mask[:,self.data.index,:]= False
		self.data.timeMask[:,self.data.index] = False
		self.data.freqMask[self.data.index,:] = False
		
		self.data.draw()
		self.data.drawSpectrum(self.data.spectrumClick)
		self.data.makeMark(self.data.spectrumClick)
		
		wx.EndBusyCursor()
		
	def onMaskResetAll(self, event):
		"""
		Reset all masks.
		"""
		
		wx.BeginBusyCursor()
		
		for i in xrange(4):
			self.data.spec.mask[:,i,:] = False
			self.data.specBandpass.mask[:,i,:]= False
			self.data.timeMask[:,i] = False
			self.data.freqMask[i,:] = False
		
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
		
		self.data.computeBandpass()
		
		if self.data.bandpass:
			self.data.draw()
			self.data.drawSpectrum(self.data.spectrumClick)
			self.data.makeMark(self.data.spectrumClick)
			
		wx.EndBusyCursor()
		
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
	
	def resizePlots(self, event):
		w, h = self.GetSize()
		dpi = self.figure1a.get_dpi()
		newW0 = 1.0*w/dpi
		newW1 = 1.0*(w/2-100)/dpi
		newW2 = 1.0*(w/2-100)/dpi
		newH0 = 1.0*(h/2-100)/dpi
		newH1 = 1.0*(h/2-200)/dpi
		self.figure1a.set_size_inches((newW1, newH0))
		self.figure1a.canvas.draw()
		self.figure1b.set_size_inches((newW2, newH0))
		self.figure1b.canvas.draw()
		self.figure2.set_size_inches((newW0, newH1))
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
		uprText = wx.TextCtrl(panel, style=wx.TE_READONLY)
		uprDec = wx.Button(panel, ID_CONTRAST_UPR_DEC, '-', size=(56, 28))
		uprInc = wx.Button(panel, ID_CONTRAST_UPR_INC, '+', size=(56, 28))
		
		lwr = wx.StaticText(panel, label='Lower Limit:')
		lwrText = wx.TextCtrl(panel, style=wx.TE_READONLY)
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
		self.bpRText.SetValue('%.2f' % self.parent.bandpassCut)
		self.dcPText.SetValue('%i'   % self.parent.driftOrder)
		self.dcCText.SetValue('%i'   % self.parent.driftCut)
		self.skSText.SetValue('%i'   % self.parent.kurtosisSec)
		self.skCText.SetValue('%i'   % self.parent.kurtosisCut)
		
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
		if self.parent.bandpassCut > 0.05:
			self.parent.bandpassCut -= 0.05
			self.bpRText.SetValue('%.2f' % self.parent.bandpassCut)
		
	def onBPRIncrease(self, event):
		if self.parent.bandpassCut < 1:
			self.parent.bandpassCut += 0.05
			self.bpRText.SetValue('%.2f' % self.parent.bandpassCut)
		
	def onDCPDecrease(self, event):
		if self.parent.driftOrder > 1:
			self.parent.driftOrder -= 1
			self.dcPText.SetValue('%i'   % self.parent.driftOrder)
		
	def onDCPIncrease(self, event):
		if self.parent.driftOrder < 12:
			self.parent.driftOrder += 1
			self.dcPText.SetValue('%i'   % self.parent.driftOrder)
	
	def onDCCDecrease(self, event):
		if self.parent.driftCut > 2:
			self.parent.driftCut -= 1
			self.dcCText.SetValue('%i'   % self.parent.driftCut)
		
	def onDCCIncrease(self, event):
		if self.parent.driftOrder < numpy.ceil(self.data.data.spec.shape[0]/300):
			self.parent.driftCut += 1
			self.dcCText.SetValue('%i'   % self.parent.driftCut)
			
	def onSKSDecrease(self, event):
		if self.parent.kurtosisSec > 1:
			self.parent.kurtosisSec -= 1
			self.skSText.SetValue('%i'   % self.parent.kurtosisSec)
		
	def onSKSIncrease(self, event):
		if self.parent.kurtosisSec < 30:
			self.parent.kurtosisSec += 1
			self.skSText.SetValue('%i'   % self.parent.kurtosisSec)
	
	def onSKCDecrease(self, event):
		if self.parent.kurtosisCut > 2:
			self.parent.kurtosisCut -= 1
			self.skCText.SetValue('%i'   % self.parent.kurtosisCut)
		
	def onSKCIncrease(self, event):
		if self.parent.kurtosisCut < 12:
			self.parent.kurtosisCut += 1
			self.skCText.SetValue('%i'   % self.parent.kurtosisCut)
	
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
		
		if self.parent.data.bandpass:
			spec = self.parent.data.specBandpass
			limits = self.parent.data.limitsBandpass
		else:
			spec = self.parent.data.spec
			limits = self.parent.data.limits
		
		# Plot Waterfall
		self.figure.clf()
		self.ax1 = self.figure.gca()
		m = self.ax1.imshow(to_dB(spec[:,self.parent.data.index,:]), interpolation='nearest', extent=(self.parent.data.freq[0]/1e6, self.parent.data.freq[-1]/1e6, self.parent.data.time[0], self.parent.data.time[-1]), origin='lower', vmin=limits[self.parent.data.index][0], vmax=limits[self.parent.data.index][1])
		cm = self.figure.colorbar(m, ax=self.ax1)
		cm.ax.set_ylabel('PSD [arb. dB]')
		self.ax1.axis('auto')
		self.ax1.set_xlabel('Frequency [MHz]')
		self.ax1.set_ylabel('Elapsed Time [s]')
		self.ax1.set_title('Tuning %i, Pol. %s' % (self.parent.data.index/2+1, 'Y' if self.parent.data.index %2 else 'X'))
		
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
			
			dataX = numpy.where(numpy.abs(clickX-self.parent.data.freq/1e6) == (numpy.abs(clickX-self.parent.data.freq/1e6).min()))[0][0]
			dataY = numpy.where(numpy.abs(clickY-self.parent.data.time) == (numpy.abs(clickY-self.parent.data.time).min()))[0][0]
			
			value = to_dB(self.parent.data.spec[dataY, self.parent.data.index, dataX])
			self.statusbar.SetStatusText("f=%.4f MHz, t=%.4f s, p=%.2f dB" % (clickX, clickY, value))
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
		
		self.drift = spec[:,:,spec.shape[2]/4:3*spec.shape[2]/4].sum(axis=2)
		
		self.ax1.plot(self.parent.data.time, self.drift[:,self.parent.data.index], linestyle='-', marker='x')
		self.ax1.set_xlim([self.parent.data.time[0], self.parent.data.time[-1]])
		self.ax1.set_xlabel('Elapsed Time [s]')
		self.ax1.set_ylabel('Total Power [arb. dB]')
		self.ax1.set_title('Tuning %i, Pol. %s' % (self.parent.data.index/2+1, 'Y' if self.parent.data.index %2 else 'X'))
		
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
			
			value = to_dB(self.drift[dataX,self.parent.data.index])
			self.statusbar.SetStatusText("t=%s, LST=%s, p=%.2f dB" % (ts, lst, value))
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
	app = wx.App(0)
	frame = MainWindow(None, -1)
	if len(args) == 1:
		frame.filename = args[0]
		frame.data = Waterfall_GUI(frame)
		frame.data.loadData(args[0])
		frame.data.draw()
	app.MainLoop()


if __name__ == "__main__":
	main(sys.argv[1:])
