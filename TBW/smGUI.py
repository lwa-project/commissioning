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

		self.color = 0
		self.antennas = antennas
		self.freq = freq
		self.spec = spec
		self.specTemplate = specTemplate
		self.resFreq = resFreq
		
		self.ax1 = None
		self.ax2 = None
		
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
		try:
			self.resFreq = dataDict['resFreq']
		except KeyError:
			self.resFreq = None

		# Set the station
		station = stations.lwa1
		self.antennas = station.getAntennas()

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
			specDiff = numpy.where( specDiff < 2, specDiff, 2)
			specDiff = numpy.where( specDiff > 0, specDiff, 0)
			
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
			specDiff = numpy.where( specDiff < 2, specDiff, 4)
			specDiff = numpy.where( specDiff > 1, specDiff, 1)
				
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
			specDiff = numpy.where( specDiff < 2, specDiff, 4)
			specDiff = numpy.where( specDiff > 1, specDiff, 1)
				
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
				specDiff = numpy.where( specDiff < 50, specDiff, 50 )
				self.resFreq = specDiff
			else:
				specDiff = self.resFreq
				specDiff = numpy.where( specDiff < 50, specDiff, 50 )
			
			cbTitle = 'Est. Resonance Point (MHz)'

		self.frame.figure.clf()
		ax1 = self.frame.figure.gca()
		# Stands 
		m = ax1.scatter(standPos[:,0], standPos[:,1]+0.8, c=specDiff[0::2], s=45.0, alpha=0.80, marker='^')
		ax1.scatter(standPos[:,0], standPos[:,1]-0.8, c=specDiff[1::2], s=45.0, alpha=0.80, marker='v')

		## Add the fence as a dashed line
		ax1.plot([-59.827, 59.771, 60.148, -59.700, -59.827], 
				[59.752, 59.864, -59.618, -59.948, 59.752], linestyle='--', color='k')

		## Add the shelter
		ax1.plot([55.863, 58.144, 58.062, 55.791, 55.863], 
				[45.946, 45.999, 51.849, 51.838, 45.946], linestyle='-', color='k')

		## Set the limits to just zoom in on the main stations
		ax1.set_xlim([-65, 65])
		ax1.set_ylim([-65, 65])

		## Set the color bar, its title, and the axis labels
		cm = self.frame.figure.colorbar(m, ax=ax1)
		cm.ax.set_ylabel(cbTitle)
		ax1.set_xlabel('$\Delta$ X [m]')
		ax1.set_ylabel('$\Delta$ Y [m]')
		
		## Draw it
		self.frame.canvas.draw()
		
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
		ax2 = self.frame.figure2.gca()
		ax2.plot(self.freq/1e6, to_dB(self.specTemplate), alpha=0.6, color='g', label='Composite')
		ax2.plot(self.freq/1e6, to_dB(self.spec[self.bestX-1,:]), color='b', label='X')
		ax2.plot(self.freq/1e6, to_dB(self.spec[self.bestY-1,:]), color='r', label='Y')
		
		## Set the axis labels and add a legend
		ax2.set_xlabel('Frequency [MHz]')
		ax2.set_ylabel('PSD [dB/RBW]')
		ax2.legend(loc=0)
		
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
		
		ax = self.frame.figure.gca()
		
		try:
			if self.oldMark is not None:
				del ax.lines[-1]
		except AttributeError:
			pass
		
		## Figure out who is who
		xy = [self.antennas[self.bestX-1].stand.x, self.antennas[self.bestX-1].stand.y]

		self.oldMark = ax.plot([xy[0], xy[0]], [xy[1], xy[1]], linestyle=' ', marker='o', ms=15.0, mfc='None', color='k')
		
		## Set the limits to just zoom in on the main stations
		ax.set_xlim([-65, 65])
		ax.set_ylim([-65, 65])


		self.frame.canvas.draw()
	
	def connect(self):
		"""
		Connect to all the events we need to interact with the plots.
		"""
		
		self.cidpress = self.frame.figure.canvas.mpl_connect('button_press_event', self.on_press)
		self.cidmotion = self.frame.figure.canvas.mpl_connect('motion_notify_event', self.on_motion)
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

			self.frame.statusbar.SetStatusText("x=%.3f m, y=%.3f m; stand #%i selected" % 
					(clickX, clickY, self.antennas[self.bestX-1].stand.id))
			
	def on_motion(self, event):
		"""
		Deal with motion events in the stand field window.  This involves 
		setting the status bar with the current x and y coordinates as well
		as the stand number of the selected stand (if any).
		"""
		
		if event.inaxes:
			clickX = event.xdata
			clickY = event.ydata
			
			if self.bestX == -1:
				# No stand selected
				self.frame.statusbar.SetStatusText("x=%.3f m, y=%.3f m" % (clickX, clickY))
			else:
				# Stand selected
				self.frame.statusbar.SetStatusText("x=%.3f m, y=%.3f m; stand #%i selected" % 
					(clickX, clickY, self.antennas[self.bestX-1].stand.id))
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
				self.frame.statusbar.SetStatusText("freq=%.3f MHz, PSD=%.3f dB/RBW; stand #%i selected" % 
					(clickX, clickY, self.antennas[self.bestX-1].stand.id))
		else:
			self.frame.statusbar.SetStatusText("")
			
	
	def disconnect(self):
		"""
		Disconnect all the stored connection ids.
		"""
		
		self.frame.figure.canvas.mpl_disconnect(self.cidpress)
		self.frame.figure.canvas.mpl_disconnect(self.cidmotion)
		self.frame.figure2.canvas.mpl_disconnect(self.cidmotion2)


ID_OPEN = 10
ID_QUIT = 11
ID_COLOR_0 = 20
ID_COLOR_1 = 21
ID_COLOR_2 = 22
ID_COLOR_3 = 23
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
		colorMenu.AppendRadioItem(ID_COLOR_0, 'Median Comparison')
		colorMenu.AppendRadioItem(ID_COLOR_3, 'Resonance Point')
		colorMenu.AppendRadioItem(ID_COLOR_1, 'RFI-46 Index')
		colorMenu.AppendRadioItem(ID_COLOR_2, 'RFI-64 Index')
		
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
		self.figure = Figure()
		self.canvas = FigureCanvasWxAgg(panel1, -1, self.figure)
		vbox1.Add(self.canvas, 1, wx.EXPAND)
		panel1.SetSizer(vbox1)
		hbox.Add(panel1, 1, wx.EXPAND)
		
		# Add a spectrum plot
		panel3 = wx.Panel(self, -1)
		vbox3 = wx.BoxSizer(wx.VERTICAL)
		self.figure2 = Figure()
		self.canvas2 = FigureCanvasWxAgg(panel3, -1, self.figure2)
		self.toolbar = NavigationToolbar2WxAgg(self.canvas2)
		self.toolbar.Realize()
		vbox3.Add(self.canvas2, 1, wx.EXPAND)
		vbox3.Add(self.toolbar, 0, wx.LEFT | wx.FIXED_MINSIZE)
		panel3.SetSizer(vbox3)
		hbox.Add(panel3, 1, wx.EXPAND)
		
		# Make the images resizable
		self.Bind(wx.EVT_PAINT, self.resizePlots)
		
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
		
		# Detail menu events
		self.Bind(wx.EVT_MENU, self.onAntenna, id=ID_DETAIL_ANT)
		self.Bind(wx.EVT_MENU, self.onStand, id=ID_DETAIL_STAND)
		self.Bind(wx.EVT_MENU, self.onFEE, id=ID_DETAIL_FEE)
		self.Bind(wx.EVT_MENU, self.onCable, id=ID_DETAIL_CABLE)
		self.Bind(wx.EVT_MENU, self.onRFI, id=ID_DETAIL_RFI)
	
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
			
	def onColor1(self, event):
		"""
		Set the stand field color coding to the RFI-46 index.  Re-draw if 
		necessary.
		"""
		
		if self.data.color != 1:
			self.data.color = 1
			self.data.draw()
			
	def onColor2(self, event):
		"""
		Set the stand field color coding to the RFI-64 index.  Re-draw if 
		necessary.
		"""
		
		if self.data.color != 2:
			self.data.color = 2
			self.data.draw()
			
	def onColor3(self, event):
		"""
		Set the stand field color coding to the estimates resonance point
		frequency in MHz.  Re-draw if necessary.
		"""
		
		if self.data.color != 3:
			self.data.color = 3
			self.data.draw()
		
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
