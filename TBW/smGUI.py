#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import numpy

from lsl.common import stations
from lsl.misc.mathutil import to_dB

import wx
import matplotlib
matplotlib.use('WXAgg')
matplotlib.interactive(True)

from matplotlib.backends.backend_wxagg import NavigationToolbar2WxAgg, FigureCanvasWxAgg
from matplotlib.figure import Figure


class TBW_GUI(object):
	"""
	Object responsible for drawing and interacting with the two matplotlib 
	canvas in the main GUI.  
	
	Arguments:
	  * frame - parent window used for displaying
	
	Keyword:
	  * antennas - list of Antenna objects that compose the stations
	  * freq - array of frequencies (in Hz) for the spectral data
	  * spec - 2-D (antennas by channels) array of spectra
	  * specTemplate - master spectrum for comparisons
	"""
	
	def __init__(self, frame, antennas=None, freq=None, spec=None, specTemplate=None, resFreq=None):
		self.frame = frame
		self.press = None
		self.bestX = -1
		self.bestY = 0

		self.filename = ''
		self.color = 0
		self.antennas = antennas
		self.freq = freq
		self.spec = spec
		self.specTemplate = specTemplate
		self.resFreq = resFreq
		
		self.ax1 = None
		self.ax2 = None
		self.oldMark = None
		
	def loadData(self, filename):
		"""
		Given the filename of an NPZ file created by stationMaster.py, load 
		in the NPZ file and set all of the various data attributes used for
		plotting.
		"""
		
		dataDict = numpy.load(filename)
		self.freq = dataDict['freq']
		masterSpectra = dataDict['masterSpectra']
		
		self.spec = masterSpectra.mean(axis=0)
		self.specTemplate = numpy.median(self.spec, axis=0)
		self.resFreq = dataDict['resFreq']

		# Set the station
		station = stations.lwa1
		self.antennas = station.getAntennas()

		# Set default colobars
		self.limits = []
		self.limits.append([0, 2])
		self.limits.append([1, 2])
		self.limits.append([1, 2])
		self.limits.append([31, 50])
		
		# Save the filename and data
		path, basename = os.path.split(filename)
		self.filename = basename
		self.date = dataDict['date']

		try:
			self.disconnect()
		except:
			pass

		# Clear the mark and plots
		self.oldMark = None
		self.frame.figure1.clf()
		self.frame.figure2.clf()

		self.connect()
	
	def draw(self):
		"""
		Draw the station stand field using the selected colorization scheme.
		The default scheme is based on the mean ratio between the spectrum 
		and the master template (self.specTemplate) between 32 and 50 MHz.
		"""
		
		wx.BeginBusyCursor()
		
		## Get Stand positions from the Antenna objects.  Select only one 
		## polarization since that is all we need
		standPos = numpy.array([[ant.stand.x, ant.stand.y, ant.stand.z] for ant in self.antennas if ant.pol == 0])

		## Get the stand quality figure-of-merit
		compLow = 32e6
		compHigh = 50e6
		
		if self.color == 0:
			# Color by mean ratio between the spectrum and the master template 
			# (self.specTemplate) between 32 and 50 MHz (default)
			specDiff = numpy.zeros(self.spec.shape[0])
			toCompare = numpy.where( (self.freq>compLow) & (self.freq<compHigh) )[0]
			for i in xrange(self.spec.shape[0]):
				specDiff[i] = (self.spec[i,toCompare] / self.specTemplate[toCompare]).mean()
			
			cbTitle = '%i to %i MHz Mean Deviation' % (compLow/1e6, compHigh/1e6)
		elif self.color == 1:
			# Color by the value of the RFI-46 index.  This index is the maximum 
			# ratio of the spectrum and the master template between 45 and 47 MHz.
			# The value of RFI-46 is also corrected for any systimatic offset 
			# between the spectrum and the template by looking at the 75 to 77 MHz
			# region.
			specDiff = numpy.zeros(self.spec.shape[0])
			rfi1 = numpy.where( (self.freq>45e6) & (self.freq<47e6) )[0]
			corr = numpy.where( (self.freq>75e6) & (self.freq<77e6) )[0]
			
			for i in xrange(self.spec.shape[0]):
				specDiff[i] = (self.spec[i,rfi1] / self.specTemplate[rfi1]).max()
				specDiff[i] /= (self.spec[i,corr] / self.specTemplate[corr]).mean()
				
			cbTitle = 'RFI-46 Index'
		elif self.color == 2:
			# Color by the value of the RFI-64 index.  This index is the maximum 
			# ratio of the spectrum and the master template between 63 and 65 MHz.
			# The value of RFI-64 is also corrected for any systimatic offset 
			# between the spectrum and the template by looking at the 75 to 77 MHz
			# region.
			specDiff = numpy.zeros(self.spec.shape[0])
			rfi2 = numpy.where( (self.freq>63e6) & (self.freq<65e6) )[0]
			corr = numpy.where( (self.freq>75e6) & (self.freq<77e6) )[0]
			
			for i in xrange(self.spec.shape[0]):
				specDiff[i] = (self.spec[i,rfi2] / self.specTemplate[rfi2]).max()
				specDiff[i] /= (self.spec[i,corr] / self.specTemplate[corr]).mean()
				
			cbTitle = 'RFI-64 Index'
		else:
			# Color by the estimated resonsonce point frequency.  This is done 
			# by finding the best-fit polynomial in between orders 3 and 12 
			# for the 31 to 70 MHz spectral region.  The best-fit polynomial is
			# then evaluated to find its maximum value and that is used as the 
			# resonance point.
			if self.resFreq is None:
				specDiff = numpy.zeros(self.spec.shape[0])
				toCompare = numpy.where( (self.freq>31e6) & (self.freq<70e6) )[0]
				for i in xrange(self.spec.shape[0]):
					bestOrder = 0
					bestRMS = 1e34
					for j in xrange(3, 12):
						coeff = numpy.polyfit(self.freq[toCompare]/1e6, to_dB(self.spec[i,toCompare]), j)
						fit = numpy.polyval(coeff, self.freq[toCompare]/1e6)
						rms = ((fit - to_dB(self.spec[i,toCompare]))**2).sum()
						if rms < bestRMS:
							bestOrder = j
							bestRMS = rms
							
					coeff = numpy.polyfit(self.freq[toCompare]/1e6, to_dB(self.spec[i,toCompare]), bestOrder)
					fit = numpy.polyval(coeff, self.freq[toCompare]/1e6)	
					specDiff[i] = self.freq[toCompare[numpy.where( fit == fit.max() )[0]]] / 1e6
				self.resFreq = specDiff
			else:
				specDiff = self.resFreq
			
			cbTitle = 'Est. Resonance Point (MHz)'

		# Clip range
		specDiff = numpy.where( specDiff < self.limits[self.color][1], specDiff, self.limits[self.color][1])
		specDiff = numpy.where( specDiff > self.limits[self.color][0], specDiff, self.limits[self.color][0])

		self.frame.figure1.clf()
		self.ax1 = self.frame.figure1.gca()
		# Stands 
		m = self.ax1.scatter(standPos[:,0], standPos[:,1]+0.8, c=specDiff[0::2], s=45.0, alpha=0.80, marker='^')
		self.ax1.scatter(standPos[:,0], standPos[:,1]-0.8, c=specDiff[1::2], s=45.0, alpha=0.80, marker='v')

		## Add the fence as a dashed line
		self.ax1.plot([-59.827, 59.771, 60.148, -59.700, -59.827], 
				[59.752, 59.864, -59.618, -59.948, 59.752], linestyle='--', color='k')

		## Add the shelter
		self.ax1.plot([55.863, 58.144, 58.062, 55.791, 55.863], 
				[45.946, 45.999, 51.849, 51.838, 45.946], linestyle='-', color='k')

		## Set the limits to just zoom in on the main station and the plot title
		if self.date is None:
			self.ax1.set_title("Filename: '%s'" % self.filename)
		else:
			self.ax1.set_title('Date: UT %s' % self.date)
		self.ax1.set_xlim([-65, 65])
		self.ax1.set_ylim([-65, 65])

		## Set the color bar, its title, and the axis labels
		cm = self.frame.figure1.colorbar(m, ax=self.ax1)
		cm.ax.set_ylabel(cbTitle)
		self.ax1.set_xlabel('$\Delta$ X [m]')
		self.ax1.set_ylabel('$\Delta$ Y [m]')
		
		if self.oldMark is not None:
			self.ax1.lines.extend(self.oldMark)
		
		## Draw it
		self.frame.canvas1.draw()
		
		wx.EndBusyCursor()

	def drawSpectrum(self, clickX, clickY):
		"""
		Get the spectra (both polarizations) for the antennas connected to 
		the selected stand.
		"""
		
		## Figure out who is who and which antennas are cloests to the 
		## clicked point.  This can be a little slow so the results are
		## saved to the bestX and bestX attributes (depending on pol.)
		dist = 1e9
		for ant in self.antennas:
			cDist = numpy.sqrt( (ant.stand.x - clickX)**2 + (ant.stand.y - clickY)**2 )
			if cDist <= dist:
				dist = cDist
				if ant.pol == 0:
					self.bestX = ant.digitizer
				else:
					self.bestY = ant.digitizer
		
		## Plot the spectra.  This plot includes the median composite 
		## (self.specTemplate) in green, the X polarization in blue, and 
		## the Y polarization in red.  
		self.frame.figure2.clf()
		self.ax2 = self.frame.figure2.gca()
		self.ax2.plot(self.freq/1e6, to_dB(self.specTemplate), alpha=0.6, color='g', label='Composite')
		self.ax2.plot(self.freq/1e6, to_dB(self.spec[self.bestX-1,:]), color='b', label='X')
		self.ax2.plot(self.freq/1e6, to_dB(self.spec[self.bestY-1,:]), color='r', label='Y')
		
		## Set the title, axis labels and add a legend
		self.ax2.set_title('Stand #%i' % self.antennas[self.bestX-1].stand.id)
		self.ax2.set_xlabel('Frequency [MHz]')
		self.ax2.set_ylabel('PSD [dB/RBW]')
		self.ax2.legend(loc=0)
		
		## Draw and save the click (Why?)
		self.frame.canvas2.draw()
		self.xClick = clickX
		self.yClick = clickY

	def makeMark(self, clickX, clickY):
		"""
		Mark the cloests stand to the clicked point.  This needs to be called
		after drawSpectrum() since that function figures out which stand is
		closest.
		"""
		
		if self.oldMark is not None:
			try:
				del self.ax1.lines[-1]
			except:
				pass
		
		## Figure out who is who
		xy = [self.antennas[self.bestX-1].stand.x, self.antennas[self.bestX-1].stand.y]

		self.oldMark = self.ax1.plot([xy[0], xy[0]], [xy[1], xy[1]], linestyle=' ', marker='o', ms=15.0, mfc='None', color='k')
		
		## Set the limits to just zoom in on the main stations
		self.ax1.set_xlim([-65, 65])
		self.ax1.set_ylim([-65, 65])


		self.frame.canvas1.draw()
	
	def connect(self):
		"""
		Connect to all the events we need to interact with the plots.
		"""
		
		self.cidpress   = self.frame.figure1.canvas.mpl_connect('button_press_event',  self.on_press)
		self.cidmotion  = self.frame.figure1.canvas.mpl_connect('motion_notify_event', self.on_motion)
		self.cidmotion2 = self.frame.figure2.canvas.mpl_connect('motion_notify_event', self.on_motion2)
	
	def on_press(self, event):
		"""
		On a button press we will update the spectrum and mark the closest 
		stand.
		"""
		
		if event.inaxes and self.frame.toolbar.mode == '':
			clickX = event.xdata
			clickY = event.ydata
				
			self.drawSpectrum(clickX, clickY)
			self.makeMark(clickX, clickY)
			
	def on_motion(self, event):
		"""
		Deal with motion events in the stand field window.  This involves 
		setting the status bar with the current x and y coordinates as well
		as the stand number of the selected stand (if any).
		"""
		
		if event.inaxes:
			clickX = event.xdata
			clickY = event.ydata
			
			self.frame.statusbar.SetStatusText("x=%.3f m, y=%.3f m" % (clickX, clickY))
		else:
			self.frame.statusbar.SetStatusText("")
			
	def on_motion2(self, event):
		"""
		Deal with motion events in the spectrum window.  This involves 
		setting the status bar with the current frequency and power values 
		as well as the stand number of the selected.  If not stand has been
		selected, nothing is shown.
		"""
		
		if event.inaxes:
			clickX = event.xdata
			clickY = event.ydata
			
			if self.bestX == -1:
				self.frame.status.bar.SetStatusText("")
			else:
				self.frame.statusbar.SetStatusText("freq=%.3f MHz, PSD=%.3f dB/RBW" % (clickX, clickY))
		else:
			self.frame.statusbar.SetStatusText("")
			
	
	def disconnect(self):
		"""
		Disconnect all the stored connection ids.
		"""
		
		self.frame.figure1.canvas.mpl_disconnect(self.cidpress)
		self.frame.figure1.canvas.mpl_disconnect(self.cidmotion)
		self.frame.figure2.canvas.mpl_disconnect(self.cidmotion2)


ID_OPEN = 10
ID_QUIT = 11
ID_COLOR_0 = 20
ID_COLOR_1 = 21
ID_COLOR_2 = 22
ID_COLOR_3 = 23
ID_COLOR_ADJUST = 24
ID_DETAIL_ANT = 30
ID_DETAIL_STAND = 31
ID_DETAIL_FEE = 32
ID_DETAIL_CABLE = 33
ID_DETAIL_RFI = 34
ID_DETAIL_SUMMARY = 35

class MainWindow(wx.Frame):
	def __init__(self, parent, id, title):
		self.dirname = ''
		self.data = None

		wx.Frame.__init__(self, parent, id, title=title, size=(600,1200))
		
		self.initUI()
		self.initEvents()
		self.Show()
		
		self.cAdjust = None
		
	def initUI(self):
		"""
		Start the user interface.
		"""
		
		self.statusbar = self.CreateStatusBar()
		
		font = wx.SystemSettings_GetFont(wx.SYS_SYSTEM_FONT)
		font.SetPointSize(10)
		
		menubar = wx.MenuBar()
		
		fileMenu = wx.Menu()
		colorMenu = wx.Menu()
		detailMenu = wx.Menu()
		
		# File menu
		open = wx.MenuItem(fileMenu, ID_OPEN, '&Open')
		fileMenu.AppendItem(open)
		fileMenu.AppendSeparator()
		quit = wx.MenuItem(fileMenu, ID_QUIT, '&Quit')
		fileMenu.AppendItem(quit)
		
		# Color menu
		colorMenu.AppendRadioItem(ID_COLOR_0, '&Median Comparison')
		colorMenu.AppendRadioItem(ID_COLOR_3, '&Resonance Point')
		colorMenu.AppendRadioItem(ID_COLOR_1, 'RFI-&46 Index')
		colorMenu.AppendRadioItem(ID_COLOR_2, 'RFI-&64 Index')
		colorMenu.AppendSeparator()
		cadj = wx.MenuItem(colorMenu, ID_COLOR_ADJUST, '&Adjust Contrast')
		colorMenu.AppendItem(cadj)
		
		# Detail menu
		dant = wx.MenuItem(detailMenu, ID_DETAIL_ANT, '&Antenna')
		detailMenu.AppendItem(dant)
		dstd = wx.MenuItem(detailMenu, ID_DETAIL_STAND, '&Stand')
		detailMenu.AppendItem(dstd)
		dfee = wx.MenuItem(detailMenu, ID_DETAIL_FEE, '&FEE')
		detailMenu.AppendItem(dfee)
		dcbl = wx.MenuItem(detailMenu, ID_DETAIL_CABLE, '&Cable')
		detailMenu.AppendItem(dcbl)
		detailMenu.AppendSeparator()
		dshl = wx.MenuItem(detailMenu, ID_DETAIL_RFI, 'Shelter &RFI Index')
		detailMenu.AppendItem(dshl)
		
		# Creating the menubar.
		menubar.Append(fileMenu, '&File')
		menubar.Append(colorMenu, '&Color Coding')
		menubar.Append(detailMenu, '&Details')
		self.SetMenuBar(menubar)
		
		hbox = wx.BoxSizer(wx.HORIZONTAL)
		
		# Add plots
		panel1 = wx.Panel(self, -1)
		vbox1 = wx.BoxSizer(wx.VERTICAL)
		self.figure1 = Figure()
		self.canvas1 = FigureCanvasWxAgg(panel1, -1, self.figure1)
		vbox1.Add(self.canvas1, 1, wx.EXPAND)
		panel1.SetSizer(vbox1)
		hbox.Add(panel1, 1, wx.EXPAND)
		
		# Add a spectrum plot
		panel2 = wx.Panel(self, -1)
		vbox2 = wx.BoxSizer(wx.VERTICAL)
		self.figure2 = Figure()
		self.canvas2 = FigureCanvasWxAgg(panel2, -1, self.figure2)
		self.toolbar = NavigationToolbar2WxAgg(self.canvas2)
		self.toolbar.Realize()
		vbox2.Add(self.canvas2, 1, wx.EXPAND)
		vbox2.Add(self.toolbar, 0, wx.LEFT | wx.FIXED_MINSIZE)
		panel2.SetSizer(vbox2)
		hbox.Add(panel2, 1, wx.EXPAND)
		
		# Use some sizers to see layout options
		self.SetSizer(hbox)
		self.SetAutoLayout(1)
		hbox.Fit(self)

	def initEvents(self):
		"""
		Set all of the various events in the main window.
		"""
		
		# File menu events
		self.Bind(wx.EVT_MENU, self.onOpen, id=ID_OPEN)
		self.Bind(wx.EVT_MENU, self.onQuit, id=ID_QUIT)
		
		# Color menu events
		self.Bind(wx.EVT_MENU, self.onColor0, id=ID_COLOR_0)
		self.Bind(wx.EVT_MENU, self.onColor1, id=ID_COLOR_1)
		self.Bind(wx.EVT_MENU, self.onColor2, id=ID_COLOR_2)
		self.Bind(wx.EVT_MENU, self.onColor3, id=ID_COLOR_3)
		self.Bind(wx.EVT_MENU, self.onAdjust, id=ID_COLOR_ADJUST)
		
		# Detail menu events
		self.Bind(wx.EVT_MENU, self.onAntenna, id=ID_DETAIL_ANT)
		self.Bind(wx.EVT_MENU, self.onStand, id=ID_DETAIL_STAND)
		self.Bind(wx.EVT_MENU, self.onFEE, id=ID_DETAIL_FEE)
		self.Bind(wx.EVT_MENU, self.onCable, id=ID_DETAIL_CABLE)
		self.Bind(wx.EVT_MENU, self.onRFI, id=ID_DETAIL_RFI)
		
		# Make the images resizable
		self.Bind(wx.EVT_PAINT, self.resizePlots)
	
	def onOpen(self,e):
		"""
		Open and load in a new NPZ file created by stationMaster.py.
		"""
		
		dlg = wx.FileDialog(self, "Choose a file", self.dirname, "", "*.*", wx.OPEN)
		if dlg.ShowModal() == wx.ID_OK:
			self.filename = dlg.GetFilename()
			self.dirname = dlg.GetDirectory()
			self.data = TBW_GUI(self)
			self.data.loadData(os.path.join(self.dirname, self.filename))
			self.data.draw()
			
			if self.cAdjust is not None:
				try:
					self.cAdjust.Close()
				except:
					pass
				self.cAdjust = None
		dlg.Destroy()
		
	def onQuit(self, event):
		"""
		Quit station master GUI.
		"""
		
		self.Close(True)
		
	def onColor0(self, event):
		"""
		Set the stand field color coding to the mean ratio of the antenna
		PSD to the median spectrum between 32 and 50 MHz.  Re-draw if 
		necessary.
		"""
		
		if self.data.color != 0:
			self.data.color = 0
			self.data.draw()
			if self.cAdjust is not None:
				try:
					self.cAdjust.Close()
				except:
					pass
			
	def onColor1(self, event):
		"""
		Set the stand field color coding to the RFI-46 index.  Re-draw if 
		necessary.
		"""
		
		if self.data.color != 1:
			self.data.color = 1
			self.data.draw()
			if self.cAdjust is not None:
				try:
					self.cAdjust.Close()
				except:
					pass
			
	def onColor2(self, event):
		"""
		Set the stand field color coding to the RFI-64 index.  Re-draw if 
		necessary.
		"""
		
		if self.data.color != 2:
			self.data.color = 2
			self.data.draw()
			if self.cAdjust is not None:
				try:
					self.cAdjust.Close()
				except:
					pass
			
	def onColor3(self, event):
		"""
		Set the stand field color coding to the estimates resonance point
		frequency in MHz.  Re-draw if necessary.
		"""
		
		if self.data.color != 3:
			self.data.color = 3
			self.data.draw()
			if self.cAdjust is not None:
				try:
					self.cAdjust.Close()
				except:
					pass
			
	def onAdjust(self, event):
		"""
		Bring up the colorbar adjustment dialog window.
		"""
		
		ContrastAdjust(self)
		
	def onAntenna(self, event):
		"""
		Display meta-data for the selected antenna pair.  This includes:
		  * ID numbers
		  * polarizations
		  * resonance point estimated from polynomial fits
		  * DP1 board number
		  * DP1 digitizer number
		  * status code
		"""
		
		if self.data is None:
			pass
		
		if self.data.bestX != -1:
			ant1 = self.data.antennas[self.data.bestX-1]
			ant2 = self.data.antennas[self.data.bestY-1]
			
			toCompare = numpy.where( (self.data.freq>31e6) & (self.data.freq<70e6) )[0]
			
			i = self.data.bestX-1
			bestOrder = 0
			bestRMS = 1e34
			for j in xrange(3, 12):
				coeff = numpy.polyfit(self.data.freq[toCompare]/1e6, to_dB(self.data.spec[i,toCompare]), j)
				fit = numpy.polyval(coeff, self.data.freq[toCompare]/1e6)
				rms = ((fit - to_dB(self.data.spec[i,toCompare]))**2).sum()
				if rms < bestRMS:
					bestOrder = j
					bestRMS = rms
						
			coeff = numpy.polyfit(self.data.freq[toCompare]/1e6, to_dB(self.data.spec[i,toCompare]), bestOrder)
			fit = numpy.polyval(coeff, self.data.freq[toCompare]/1e6)	
			res1 = self.data.freq[toCompare[numpy.where( fit == fit.max() )[0]]] / 1e6
			
			i = self.data.bestY-1
			bestOrder = 0
			bestRMS = 1e34
			for j in xrange(3, 12):
				coeff = numpy.polyfit(self.data.freq[toCompare]/1e6, to_dB(self.data.spec[i,toCompare]), j)
				fit = numpy.polyval(coeff, self.data.freq[toCompare]/1e6)
				rms = ((fit - to_dB(self.data.spec[i,toCompare]))**2).sum()
				if rms < bestRMS:
					bestOrder = j
					bestRMS = rms
						
			coeff = numpy.polyfit(self.data.freq[toCompare]/1e6, to_dB(self.data.spec[i,toCompare]), bestOrder)
			fit = numpy.polyval(coeff, self.data.freq[toCompare]/1e6)	
			res2 = self.data.freq[toCompare[numpy.where( fit == fit.max() )[0]]] / 1e6
			
			outString = """Antenna: %i
Polarization: %i

Est. Resonance: %.3f MHz

DP1 Board: %i
Digitizer: %i

Status: %i

---

Antenna: %i
Polarization: %i

Est. Resonance: %.3f MHz

DP1 Board: %i
Digitizer: %i

Status: %i
""" % (ant1.id, ant1.pol, res1, ant1.board, ant1.digitizer, ant1.status, 
		ant2.id, ant2.pol, res2, ant2.board, ant2.digitizer, ant2.status)
		
			box = wx.MessageDialog(self, outString, "Antenna Details")
			box.ShowModal()
		else:
			pass
	
	def onStand(self, event):
		"""
		Display meta-data about the selected stand.  This includes:
		  * ID number
		  * position
		  * distance from the center of the shelter
		  * distance from the fence
		"""
		
		if self.data is None:
			pass
		
		if self.data.bestX != -1:
			std = self.data.antennas[self.data.bestX-1].stand
			shlDist = numpy.sqrt( (std.x - 56.965)**2 + (std.y - 48.908)**2 )
			
			fenDistA = numpy.zeros(4)
			k = 0
			for p1,p2 in zip([(-59.827,59.752), (59.771,59.864), (60.148,-59.618), (-59.700,-59.948)], [(59.771,59.864), (60.148,-59.618), (-59.700,-59.948), (-59.827,59.752), (59.771,59.864)]):
				x1 = p1[0]
				y1 = p1[1]
				x2 = p2[0]
				y2 = p2[1]
				
				a = (y2-y1)/(x2-x1)
				b = (y2*x1-y1*x2)/(x2-x1)
				
				x3 = std.x
				y3 = std.y
				
				x4 = (x3/a + y3 - b)*a / (a**2+1)
				y4 = a*x4 + b
				
				fenDistA[k] = numpy.sqrt( (x3-x4)**2 + (y3-y4)**2 )
				k += 1
				
			fenDist = fenDistA.min()
			
			outString = """Stand: %i

Coordinates:
x = %.3f m
y = %.3f m
z = %.3f m

Shelter Distance: %.3f m
Fence Distance: %.3f m
""" % (std.id, std.x, std.y, std.z, shlDist, fenDist)
			
			box = wx.MessageDialog(self, outString, "Stand Detail")
			box.ShowModal()
		else:
			pass
	
	def onCable(self, event):
		"""
		Display meta-data about the RPD cables connecting the the selected 
		stand/antennas back to the shelter.  This includes:
		  * ID names
		  * lengths
		  * delay at 10 MHz
		  * gain at 10 MHz
		"""
		
		if self.data is None:
			pass
		
		if self.data.bestX != -1:
			ant1 = self.data.antennas[self.data.bestX-1]
			ant2 = self.data.antennas[self.data.bestY-1]
			
			outString = """Antenna: %i
Cable: %s

Length: %.1f m

Delay @ 10 MHz: %.1f ns
Gain @ 10 MHz: %.1f dB

---

Antenna: %i
Cable: %s

Length: %.1f m

Delay @ 10 MHz: %.1f ns
Gain @ 10 MHz: %.1f dB
""" % (ant1.id, ant1.cable.id, ant1.cable.length, ant1.cable.delay(10e6, ns=True), to_dB(ant1.cable.gain(10e6)), 
		ant2.id, ant2.cable.id, ant2.cable.length, ant2.cable.delay(10e6, ns=True), to_dB(ant2.cable.gain(10e6)))
			
			box = wx.MessageDialog(self, outString, "Cable Details")
			box.ShowModal()
		else:
			pass
		
	def onFEE(self, event):
		"""
		Display meta-data about the FEE installed in the selecte stand.  
		This includes:
		  * ID name
		  * which antennas are connected to which ports
		  * gain settings
		  * status
		"""
		
		if self.data is None:
			pass
		
		if self.data.bestX != -1:
			fee = self.data.antennas[self.data.bestX-1].fee
			portXA = self.data.antennas[self.data.bestX-1].id
			portXP = self.data.antennas[self.data.bestX-1].feePort
			portYA = self.data.antennas[self.data.bestY-1].id
			portYP = self.data.antennas[self.data.bestY-1].feePort
			
			outString = """FEE: %s

Ports:
%i = antenna %i (X pol.)
%i = antenna %i (Y pol.)

Gains:
1 = %.3f dB
2 = %.3f dB

Status: %i
""" % (fee.id, portXP, portXA, portYP, portYA, fee.gain1, fee.gain2, fee.status)
			
			box = wx.MessageDialog(self, outString, "FEE Detail")
			box.ShowModal()
		else:
			pass
		
	def onRFI(self, event):
		if self.data is None:
			pass
		
		if self.data.bestX != -1:
			ant1 = self.data.antennas[self.data.bestX-1]
			ant2 = self.data.antennas[self.data.bestY-1]
			
			rfi1 = numpy.where( (self.data.freq>45e6) & (self.data.freq<47e6) )[0]
			rfi2 = numpy.where( (self.data.freq>63e6) & (self.data.freq<65e6) )[0]
			corr = numpy.where( (self.data.freq>75e6) & (self.data.freq<77e6) )[0]
			
			a1r1 = (self.data.spec[self.data.bestX-1,rfi1] / self.data.specTemplate[rfi1]).max()
			a1r2 = (self.data.spec[self.data.bestX-1,rfi2] / self.data.specTemplate[rfi2]).max()
			c1 = (self.data.spec[self.data.bestX-1,corr] / self.data.specTemplate[corr]).mean()
			a2r1 = (self.data.spec[self.data.bestY-1,rfi1] / self.data.specTemplate[rfi1]).max()
			a2r2 = (self.data.spec[self.data.bestY-1,rfi2] / self.data.specTemplate[rfi2]).max()
			c2 = (self.data.spec[self.data.bestY-1,corr] / self.data.specTemplate[corr]).mean()
			
			outString = """Antenna: %i
Polarization: %i

RFI-46 Index:
raw = %.3f
corrected = %.3f

RFI-64 Index:
raw = %.3f
corrected = %.3f

---

Antenna: %i
Polarization: %i

RFI-46 Index:
raw = %.3f
corrected = %.3f

RFI-64 Index:
raw = %.3f
corrected = %.3f
""" % (ant1.id, ant1.pol, a1r1, a1r1/c1, a1r2, a1r2/c1, 
		ant2.id, ant2.pol, a2r1, a2r1/c2, a2r2, a2r2/c2)
		
			box = wx.MessageDialog(self, outString, "Shelter RFI Details")
			box.ShowModal()
		else:
			pass

	def resizePlots(self, event):
		w, h = self.GetSize()
		dpi = self.figure.get_dpi()
		newW = 1.0*w/dpi
		newH1 = 1.0*(h/2-100)/dpi
		newH2 = 1.0*(h/2-75)/dpi
		self.figure.set_size_inches((newW, newH1))
		self.figure.canvas.draw()
		self.figure2.set_size_inches((newW, newH2))
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
		
		if self.parent.data.color == 0:
			mode = 'Median Comparision'
		elif self.parent.data.color == 1:
			mode = 'RFI-46 Index'
		elif self.parent.data.color == 2:
			mode = 'RFI-64 Index'
		else:
			mode = 'Resonance Point'
		typ = wx.StaticText(panel, label='Color Coding Mode: %s' % mode)
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
		color = self.parent.data.color
		self.uText.SetValue('%.1f' % self.parent.data.limits[color][1])
		self.lText.SetValue('%.1f' % self.parent.data.limits[color][0])
		self.rText.SetValue('%.1f' % self.__getRange(color))
		
	def initEvents(self):
		self.Bind(wx.EVT_BUTTON, self.onUpperDecrease, id=ID_CONTRAST_UPR_DEC)
		self.Bind(wx.EVT_BUTTON, self.onUpperIncrease, id=ID_CONTRAST_UPR_INC)
		self.Bind(wx.EVT_BUTTON, self.onLowerDecrease, id=ID_CONTRAST_LWR_DEC)
		self.Bind(wx.EVT_BUTTON, self.onLowerIncrease, id=ID_CONTRAST_LWR_INC)
		
		self.Bind(wx.EVT_BUTTON, self.onOk, id=ID_CONTRAST_OK)
		
	def onUpperDecrease(self, event):
		color = self.parent.data.color
		self.parent.data.limits[color][1] -= self.__getIncrement(color)
		self.uText.SetValue('%.1f' % self.parent.data.limits[color][1])
		self.rText.SetValue('%.1f' % self.__getRange(color))
		self.parent.data.draw()
		
	def onUpperIncrease(self, event):
		color = self.parent.data.color
		self.parent.data.limits[color][1] += self.__getIncrement(color)
		self.uText.SetValue('%.1f' % self.parent.data.limits[color][1])
		self.rText.SetValue('%.1f' % self.__getRange(color))
		self.parent.data.draw()
		
	def onLowerDecrease(self, event):
		color = self.parent.data.color
		self.parent.data.limits[color][0] -= self.__getIncrement(color)
		self.lText.SetValue('%.1f' % self.parent.data.limits[color][0])
		self.rText.SetValue('%.1f' % self.__getRange(color))
		self.parent.data.draw()
		
	def onLowerIncrease(self, event):
		color = self.parent.data.color
		self.parent.data.limits[color][0] += self.__getIncrement(color)
		self.lText.SetValue('%.1f' % self.parent.data.limits[color][0])
		self.rText.SetValue('%.1f' % self.__getRange(color))
		self.parent.data.draw()
		
	def onOk(self, event):
		self.parent.cAdjust = None
		self.Close()
	
	def __getRange(self, color):
		return (self.parent.data.limits[color][1] - self.parent.data.limits[color][0])
		
	def __getIncrement(self, color):
		return 0.1*self.__getRange(color)
	

def main(args):
	app = wx.App(0)
	frame = MainWindow(None, -1, "Station Master GUI")
	if len(args) == 1:
		frame.data = TBW_GUI(frame)
		frame.data.loadData(args[0])
		frame.data.draw()
	app.MainLoop()


if __name__ == "__main__":
	main(sys.argv[1:])
