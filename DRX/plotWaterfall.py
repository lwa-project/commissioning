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
import numpy
from datetime import datetime

from lsl.common.dp import fS
from lsl.misc.mathutil import to_dB

import wx
import matplotlib
matplotlib.use('WXAgg')
matplotlib.interactive(True)

from matplotlib.backends.backend_wxagg import NavigationToolbar2WxAgg, FigureCanvasWxAgg
from matplotlib.figure import Figure

try:
	from dp import drxFilter
except ImportError:
	from scipy.signal import freqz
	from scipy.interpolate import interp1d
	
	_nPts = 500
	
	def drxFilter(sampleRate=19.6e6):
		"""
		Return a function that will generate the shape of a DRX filter for a given sample
		rate.
		
		Based on memo DRX0001.
		"""
		
		decimation = fS / sampleRate
		decimationCIC = decimation / 2
		
		# CIC settings
		N = 5
		R = 5
		
		# FIR coefficients
		drxFIR = [-6.2000000000000000e+001,  6.6000000000000000e+001,  1.4500000000000000e+002, 
				3.4000000000000000e+001, -1.4400000000000000e+002, -5.9000000000000000e+001, 
				1.9900000000000000e+002,  1.4500000000000000e+002, -2.2700000000000000e+002, 
				-2.5700000000000000e+002,  2.3200000000000000e+002,  4.0500000000000000e+002, 
				-1.9400000000000000e+002, -5.8300000000000000e+002,  9.2000000000000000e+001, 
				7.8200000000000000e+002,  9.4000000000000000e+001, -9.9000000000000000e+002, 
				-3.9700000000000000e+002,  1.1860000000000000e+003,  8.5900000000000000e+002, 
				-1.3400000000000000e+003, -1.5650000000000000e+003,  1.3960000000000000e+003, 
				2.7180000000000000e+003, -1.1870000000000000e+003, -4.9600000000000000e+003, 
				-1.8900000000000000e+002,  1.1431000000000000e+004,  1.7747000000000000e+004, 
				1.1431000000000000e+004, -1.8900000000000000e+002, -4.9600000000000000e+003, 
				-1.1870000000000000e+003,  2.7180000000000000e+003,  1.3960000000000000e+003, 
				-1.5650000000000000e+003, -1.3400000000000000e+003,  8.5900000000000000e+002, 
				1.1860000000000000e+003, -3.9700000000000000e+002, -9.9000000000000000e+002,
				9.4000000000000000e+001,  7.8200000000000000e+002,  9.2000000000000000e+001,
				-5.8300000000000000e+002, -1.9400000000000000e+002,  4.0500000000000000e+002, 
				2.3200000000000000e+002, -2.5700000000000000e+002, -2.2700000000000000e+002, 
				1.4500000000000000e+002,  1.9900000000000000e+002, -5.9000000000000000e+001, 
				-1.4400000000000000e+002,  3.4000000000000000e+001,  1.4500000000000000e+002, 
				6.6000000000000000e+001, -6.2000000000000000e+001]
			
		# Part 1 - CIC filter
		h = numpy.linspace(0, numpy.pi/decimationCIC/2, num=_nPts, endpoint=True)
		wCIC = (numpy.sin(h*R)/numpy.sin(h/2))**N
		wCIC[0] = (2*R)**N
		
		# Part 2 - FIR filter
		h, wFIR = freqz(drxFIR, 1, _nPts)
		
		# Cascade
		w = numpy.abs(wCIC) * numpy.abs(wFIR)
		
		# Convert to a "real" frequency and magnitude response
		h *= fS / decimation / numpy.pi
		w = numpy.abs(w)**2
		
		# Mirror
		h = numpy.concatenate([-h[::-1], h[1:]])
		w = numpy.concatenate([w[::-1], w[1:]])
		
		# Return the interpolating function
		return interp1d(h, w/w.max(), kind='cubic', bounds_error=False)



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
		
		# Save the filename
		self.filename = filename
		dataDictionary = numpy.load(filename)
		
		# Load the Data
		try:
			self.srate = dataDictionary['srate']
		except KeyError:
			self.srate = 19.6e6
		self.freq = dataDictionary['freq']
		try:
			mask = dataDictionary['mask']
		except KeyError:
			mask = numpy.zeros(dataDictionary['spec'].shape, dtype=numpy.int16)
		self.spec = numpy.ma.array(dataDictionary['spec'], mask=mask)
		self.tInt = dataDictionary['tInt']
		self.time = self.tInt * numpy.arange(self.spec.shape[0])
		
		# Other data to keep around in case we save
		self.timesNPZ = dataDictionary['times']
		self.standMapperNPZ = dataDictionary['standMapper']
		
		# Get the filter model
		self.bpm = drxFilter(sampleRate=self.srate)(self.freq)
		
		# Compute the bandpass fit
		self.computeBandpass()
		
		# Set default colobars
		self.limits = []
		for i in xrange(self.spec.shape[1]):
			self.limits.append( [to_dB(self.spec[:,i,:]).min(), to_dB(self.spec[:,i,:]).max()] )
			
		self.limitsBandpass = []
		toUse = numpy.arange(self.spec.shape[2]/10, 9*self.spec.shape[2]/10)
		for i in xrange(self.spec.shape[1]):
			self.limitsBandpass.append( [to_dB(self.specBandpass[:,i,toUse]).min(), to_dB(self.specBandpass[:,i,toUse]).max()] )
		
		try:
			self.diconnect()
		except:
			pass
		
		# Clear the old marks
		self.oldMarkA = None
		self.oldMarkB = None
		self.frame.figure1a.clf()
		self.frame.figure1a.clf()
		self.frame.figure2.clf()
		
		self.connect()
	
	def computeBandpass(self):
		"""
		Compute the bandpass fits.
		"""
		
		# Account for the ARX bandpass by fitting to the inner 80% of the band
		toUse = numpy.arange(self.spec.shape[2]/10, 9*self.spec.shape[2]/10)
		
		medianSpec = numpy.median(self.spec, axis=0)
		bpm2 = []
		for i in xrange(self.spec.shape[1]):
			coeff = numpy.polyfit(self.freq[toUse], self.bpm[toUse]/self.bpm.mean() * medianSpec[i,toUse].mean()/medianSpec[i,toUse], 9)
			noiseSlope = numpy.polyval(coeff, self.freq)
			bpm2.append( self.bpm / noiseSlope )
			
		# Apply the bandpass correction
		self.specBandpass = numpy.ma.array(self.spec.data*1.0, mask=self.spec.mask)
		for i in xrange(self.spec.shape[1]):
			for j in xrange(self.spec.shape[0]):
				self.specBandpass[j,i,:] = bpm2[i].mean() / bpm2[i] * self.spec[j,i,:]
	
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
		m = self.ax1a.imshow(to_dB(spec[:,self.index,:]), extent=(self.freq[0]/1e6, self.freq[-1]/1e6, self.time[0], self.time[-1]), origin='lower', vmin=limits[self.index][0], vmax=limits[self.index][1])
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
		self.drift = spec[:,:,256:768].sum(axis=2)
		
		self.frame.figure1b.clf()
		self.ax1b = self.frame.figure1b.gca()
		self.ax1b.plot(self.drift[:,self.index], self.time, linestyle=' ', marker='x')
		self.ax1b.set_ylim([self.time[0], self.time[-1]])
		self.ax1b.set_xlabel('Total Power [arb. dB]')
		self.ax1b.set_ylabel('Elapsed Time [s]')
		
		if self.oldMarkB is not None:
			self.ax1b.lines.extend(self.oldMarkB)
		
		self.frame.canvas1b.draw()
	
	def drawSpectrum(self, clickY):
		"""Get the spectrum at a particular point in time."""
		
		if self.bandpass:
			spec = self.specBandpass[:,self.index,:]
			limits = self.limitsBandpass
		else:
			spec = self.spec[:,self.index,:]
			limits = self.limits
		medianSpec = numpy.median(spec, axis=0)
		
		if self.frame.toolbar.mode == 'zoom rect':
			try:
				oldXlim = self.ax2.get_xlim()
				oldYlim = self.ax2.get_ylim()
			except:
				oldXlim = [self.freq[0]/1e6, self.freq[-1]/1e6]
				oldYlim = limits[self.index]
		else:
			oldXlim = [self.freq[0]/1e6, self.freq[-1]/1e6]
			oldYlim = limits[self.index]
			
		dataY = numpy.where(numpy.abs(clickY-self.time) == (numpy.abs(clickY-self.time).min()))[0][0]
		
		self.frame.figure2.clf()
		self.ax2 = self.frame.figure2.gca()
		self.ax2.plot(self.freq/1e6, spec[dataY,:], linestyle=' ', marker='o', label='Current', color='blue')
		self.ax2.plot(self.freq/1e6, medianSpec, label='Median', alpha=0.5, color='green')
		self.ax2.set_xlim(oldXlim)
		self.ax2.set_ylim(oldYlim)
		self.ax2.legend(loc=0)
		self.ax2.set_xlabel('Frequency [MHz]')
		self.ax2.set_ylabel('PSD [arb. dB]')
		
		if self.bandpass:
			self.ax2.set_title("%s UTC + bandpass" % datetime.utcfromtimestamp(self.timesNPZ[dataY]))
		else:
			self.ax2.set_title("%s UTC" % datetime.utcfromtimestamp(self.timesNPZ[dataY]))
		
		self.frame.canvas2.draw()
		self.spectrumClick = clickY
	
	def makeMark(self, clickY):
		
		dataY = numpy.where(numpy.abs(clickY-self.time) == (numpy.abs(clickY-self.time).min()))[0][0]
		
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
	
	def on_press1a(self, event):
		"""
		On button press we will see if the mouse is over us and store some data
		"""
		
		if event.inaxes:
			clickX = event.xdata
			clickY = event.ydata
			
			dataY = numpy.where(numpy.abs(clickY-self.time) == (numpy.abs(clickY-self.time).min()))[0][0]
			
			if event.button == 1:
				self.drawSpectrum(clickY)
				self.makeMark(clickY)
			elif event.button == 3:
				self.spec.mask[dataY, self.index, :] = 1
				self.specBandpass.mask[dataY, self.index, :] = 1
				self.draw()
			elif event.button == 2:
				self.spec.mask[dataY, self.index, :] = 0
				self.specBandpass.mask[dataY, self.index, :] = 0
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
			
			d = numpy.sqrt( ((clickX - self.drift.data[:,self.index])/rangeX)**2 + ((clickY - self.time)/rangeY)**2 )
			best = numpy.where( d == d.min() )
			bestD = d[best]
					
			print best, bestD, clickX, clickY, self.drift.data[best,self.index], self.time[best]
			
			if event.button == 1:
				self.drawSpectrum(clickY)
				self.makeMark(clickY)
			elif event.button == 3:
				self.spec.mask[best, self.index, :] = 1
				self.specBandpass.mask[best, self.index, :] = 1
				self.draw()
			elif event.button == 2:
				self.spec.mask[best, self.index, :] = 0
				self.specBandpass.mask[best, self.index, :] = 0
				self.draw()
			else:
				pass
			
	def on_press2(self, event):
		if event.inaxes:
			clickX = event.xdata
			clickY = event.ydata
			
			dataX = numpy.where(numpy.abs(clickX-self.freq/1e6) == (numpy.abs(clickX-self.freq/1e6).min()))[0][0]
			
			if event.button == 3:
				self.spec.mask[:, self.index, dataX] = 1
				self.specBandpass.mask[:, self.index, dataX] = 1
			elif event.button == 2:
				self.spec.mask[:, self.index, dataX] = 0
				self.specBandpass.mask[:, self.index, dataX] = 0
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
ID_COLOR_ADJUST = 20
ID_TUNING1_X = 30
ID_TUNING1_Y = 31
ID_TUNING2_X = 32
ID_TUNING2_Y = 33
ID_BANDPASS_ON = 40
ID_BANDPASS_OFF = 41
ID_BANDPASS_RECOMPUTE = 42

class MainWindow(wx.Frame):
	def __init__(self, parent, id):
		self.dirname = ''
		self.filename = ''
		self.data = None
		
		wx.Frame.__init__(self, parent, id, title="DRX Waterfall Viewer", size=(600,800))
		
		self.initUI()
		self.initEvents()
		self.Show()
		
		self.cAdjust = None
		
	def initUI(self):
		self.statusbar = self.CreateStatusBar() # A Statusbar in the bottom of the window
		
		font = wx.SystemSettings_GetFont(wx.SYS_SYSTEM_FONT)
		font.SetPointSize(10)
		
		menuBar = wx.MenuBar()
		
		fileMenu = wx.Menu()
		colorMenu = wx.Menu()
		dataMenu = wx.Menu()
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
		cadj = wx.MenuItem(colorMenu, ID_COLOR_ADJUST, '&Adjust Contrast')
		colorMenu.AppendItem(cadj)
		
		## Data Menu
		dataMenu.AppendRadioItem(ID_TUNING1_X, 'Tuning 1, Pol. X')
		dataMenu.AppendRadioItem(ID_TUNING1_Y, 'Tuning 1, Pol. Y')
		dataMenu.AppendSeparator()
		dataMenu.AppendRadioItem(ID_TUNING2_X, 'Tuning 2, Pol. X')
		dataMenu.AppendRadioItem(ID_TUNING2_Y, 'Tuning 2, Pol. Y')
		
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
		
		self.Bind(wx.EVT_MENU, self.onAdjust, id=ID_COLOR_ADJUST)
		
		self.Bind(wx.EVT_MENU, self.onTuning1X, id=ID_TUNING1_X)
		self.Bind(wx.EVT_MENU, self.onTuning1Y, id=ID_TUNING1_Y)
		self.Bind(wx.EVT_MENU, self.onTuning2X, id=ID_TUNING2_X)
		self.Bind(wx.EVT_MENU, self.onTuning2Y, id=ID_TUNING2_Y)
		
		self.Bind(wx.EVT_MENU, self.onBandpassOn, id=ID_BANDPASS_ON)
		self.Bind(wx.EVT_MENU, self.onBandpassOff, id=ID_BANDPASS_OFF)
		self.Bind(wx.EVT_MENU, self.onBandpassRecompute, id=ID_BANDPASS_RECOMPUTE)
		
		# Key events
		self.canvas1a.Bind(wx.EVT_KEY_UP, self.onKeyPress)
		self.canvas1b.Bind(wx.EVT_KEY_UP, self.onKeyPress)
		
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
			if self.data.spectrumClick is not None:
				self.data.drawSpectrum(self.data.spectrumClick + self.data.tInt)
				self.data.makeMark(self.data.spectrumClick + self.data.tInt)
		elif keycode == wx.WXK_DOWN:
			if self.data.spectrumClick is not None:
				self.data.drawSpectrum(self.data.spectrumClick - self.data.tInt)
				self.data.makeMark(self.data.spectrumClick - self.data.tInt)
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
		self.parent.data.limits[index][1] -= self.__getIncrement(index)
		self.uText.SetValue('%.1f' % self.parent.data.limits[index][1])
		self.rText.SetValue('%.1f' % self.__getRange(index))
		self.parent.data.draw()
		
	def onUpperIncrease(self, event):
		index = self.parent.data.index
		self.parent.data.limits[index][1] += self.__getIncrement(index)
		self.uText.SetValue('%.1f' % self.parent.data.limits[index][1])
		self.rText.SetValue('%.1f' % self.__getRange(index))
		self.parent.data.draw()
		
	def onLowerDecrease(self, event):
		index = self.parent.data.index
		self.parent.data.limits[index][0] -= self.__getIncrement(index)
		self.lText.SetValue('%.1f' % self.parent.data.limits[index][0])
		self.rText.SetValue('%.1f' % self.__getRange(index))
		self.parent.data.draw()
		
	def onLowerIncrease(self, event):
		index = self.parent.data.index
		self.parent.data.limits[index][0] += self.__getIncrement(index)
		self.lText.SetValue('%.1f' % self.parent.data.limits[index][0])
		self.rText.SetValue('%.1f' % self.__getRange(index))
		self.parent.data.draw()
		
	def onOk(self, event):
		self.parent.cAdjust = None
		self.Close()
	
	def __getRange(self, index):
		return (self.parent.data.limits[index][1] - self.parent.data.limits[index][0])
		
	def __getIncrement(self, index):
		return 0.1*self.__getRange(index)


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
