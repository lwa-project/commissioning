#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Display NPZ data from stationMaster in an interactive GUI sort of way.

$Rev$
$LastChangedBy$
$LastChangedDate$
"""

import os
import sys
import numpy
import tempfile

from lsl.common import stations
from lsl.misc.mathutil import to_dB

import wx
import matplotlib
matplotlib.use('WXAgg')
matplotlib.interactive(True)

from matplotlib.backends.backend_wxagg import NavigationToolbar2WxAgg, FigureCanvasWxAgg
from matplotlib.figure import Figure

class PlotPanel(wx.Panel):
	"""
	The PlotPanel has a Figure and a Canvas. OnSize events simply set a 
	flag, and the actual resizing of the figure is triggered by an Idle event.
	
	From: http://www.scipy.org/Matplotlib_figure_in_a_wx_panel
	"""
	
	def __init__(self, parent, color=None, dpi=None, **kwargs):
		from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg
		from matplotlib.figure import Figure

		# initialize Panel
		if 'id' not in kwargs.keys():
			kwargs['id'] = wx.ID_ANY
		if 'style' not in kwargs.keys():
			kwargs['style'] = wx.NO_FULL_REPAINT_ON_RESIZE
		wx.Panel.__init__(self, parent, **kwargs)

		# initialize matplotlib stuff
		self.figure = Figure(None, dpi)
		self.canvas = FigureCanvasWxAgg(self, -1, self.figure)
		self.SetColor(color)

		self._SetSize()
		self.draw()

		self._resizeflag = False

		self.Bind(wx.EVT_IDLE, self._onIdle)
		self.Bind(wx.EVT_SIZE, self._onSize)

	def SetColor( self, rgbtuple=None ):
		"""
		Set figure and canvas colours to be the same.
		"""
		
		if rgbtuple is None:
			rgbtuple = wx.SystemSettings.GetColour( wx.SYS_COLOUR_BTNFACE ).Get()
		clr = [c/255. for c in rgbtuple]
		self.figure.set_facecolor(clr)
		self.figure.set_edgecolor(clr)
		self.canvas.SetBackgroundColour(wx.Colour(*rgbtuple))

	def _onSize(self, event):
		self._resizeflag = True

	def _onIdle(self, evt):
		if self._resizeflag:
			self._resizeflag = False
			self._SetSize()

	def _SetSize(self):
		pixels = tuple(self.parent.GetClientSize())
		self.SetSize(pixels)
		self.canvas.SetSize(pixels)
		self.figure.set_size_inches(float( pixels[0] )/self.figure.get_dpi(), float( pixels[1] )/self.figure.get_dpi())

	def draw(self):
		pass # abstract, to be overridden by child classes


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
		self.avgPower = None
		self.dataRange = None
		
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
		try:
			self.avgPower = dataDict['avgPower']
		except KeyError:
			self.avgPower = None
		try:
			self.dataRange = dataDict['dataRange']
		except KeyError:
			self.dataRange = None

		# Set the station
		try:
			ssmifContents = dataDict['ssmifContents']
			if ssmifContents.shape == ():
				station = stations.lwa2
				antennas = station.getAntennas()
			else:
				fh, tempSSMIF = tempfile.mkstemp(suffix='.txt', prefix='ssmif-')
				fh = open(tempSSMIF, 'w')
				for line in ssmifContents:
					fh.write('%s\n' % line)
				fh.close()
				
				station = stations.parseSSMIF(tempSSMIF)
				antennas = station.getAntennas()
				os.unlink(tempSSMIF)
			
		except KeyError:
			station = stations.lwa2
			antennas = station.getAntennas()
		self.antennas = []
		for a in antennas:
			if a.digitizer != 0:
				self.antennas.append(a)
		del(antennas)

		# Set default colobars
		self.limits = []
		self.limits.append([0, 2])
		self.limits.append([1, 2])
		self.limits.append([1, 2])
		self.limits.append([0, 3])
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
		elif self.color == 3:
			# Color by antenna status code.
			specDiff = numpy.zeros(self.spec.shape[0])
			for i in xrange(self.spec.shape[0]):
				specDiff[i] = self.antennas[i].status
				
			cbTitle = 'Antenna Status'
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
		self.ax1.plot([-44.315, 72.150, 44.077, -72.543, -44.315], 
				[-72.522, -44.277, 72.191, 43.972, -72.522], linestyle='--', color='k')

		### Add the shelter
		#self.ax1.plot([55.863, 58.144, 58.062, 55.791, 55.863], 
				#[45.946, 45.999, 51.849, 51.838, 45.946], linestyle='-', color='k')

		## Set the limits to just zoom in on the main station and the plot title
		if self.date is None:
			self.ax1.set_title("Filename: '%s'" % self.filename)
		else:
			self.ax1.set_title('Date: UT %s' % self.date)
		self.ax1.set_xlim([-75, 75])
		self.ax1.set_ylim([-75, 75])

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

	def drawSpectrum(self, clickX, clickY, preferStand=None):
		"""
		Get the spectra (both polarizations) for the antennas connected to 
		the selected stand.
		"""
		
		if preferStand is None:
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
		
		else:
			## Right now 259 and 260 are at 0,0,0 and sit on top of each other.  Using
			## the preferStand keyword, we can break this at least for searches
			for ant in self.antennas:
				if preferStand == ant.stand.id:
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
		self.ax1.set_xlim([-75, 75])
		self.ax1.set_ylim([-75, 75])


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
ID_SSMIF = 11
ID_QUIT = 12
ID_COLOR_0 = 20
ID_COLOR_1 = 21
ID_COLOR_2 = 22
ID_COLOR_3 = 23
ID_COLOR_4 = 24
ID_COLOR_ADJUST = 25
ID_DETAIL_ANT = 30
ID_DETAIL_STAND = 31
ID_DETAIL_FEE = 32
ID_DETAIL_CABLE = 33
ID_DETAIL_RFI = 34
ID_DETAIL_CHANGE_STATUS = 35
ID_AVG_POWER = 40
ID_AVG_RANGE = 41
ID_AVG_SUMMARY = 42
ID_SELECT_DIGITIZER = 50
ID_SELECT_ANTENNA = 51
ID_SELECT_STAND = 52

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
		powerMenu = wx.Menu()
		selectMenu = wx.Menu()
		
		## File menu
		open = wx.MenuItem(fileMenu, ID_OPEN, '&Open')
		fileMenu.AppendItem(open)
		ssmif = wx.MenuItem(fileMenu, ID_SSMIF, '&Show SSMIF Status')
		fileMenu.AppendItem(ssmif)
		fileMenu.AppendSeparator()
		quit = wx.MenuItem(fileMenu, ID_QUIT, '&Quit')
		fileMenu.AppendItem(quit)
		
		# Color menu
		colorMenu.AppendRadioItem(ID_COLOR_0, '&Median Comparison')
		colorMenu.AppendRadioItem(ID_COLOR_4, '&Resonance Point')
		colorMenu.AppendRadioItem(ID_COLOR_1, 'RFI-&46 Index')
		colorMenu.AppendRadioItem(ID_COLOR_2, 'RFI-&64 Index')
		colorMenu.AppendRadioItem(ID_COLOR_3, 'Antenna Status')
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
		detailMenu.AppendSeparator()
		dcst = wx.MenuItem(detailMenu, ID_DETAIL_CHANGE_STATUS, 'Change Antenna/FEE Status')
		detailMenu.AppendItem(dcst)
		
		# Power
		apwr = wx.MenuItem(powerMenu, ID_AVG_POWER, '&Plot Power')
		powerMenu.AppendItem(apwr)
		drng = wx.MenuItem(powerMenu, ID_AVG_RANGE, '&Plot Data Range')
		powerMenu.AppendItem(drng)
		powerMenu.AppendSeparator()
		spwr = wx.MenuItem(powerMenu, ID_AVG_SUMMARY, '&Summary')
		powerMenu.AppendItem(spwr)
		
		# Select
		sant = wx.MenuItem(selectMenu, ID_SELECT_ANTENNA, '&Antenna ID')
		selectMenu.AppendItem(sant)
		sstd = wx.MenuItem(selectMenu, ID_SELECT_STAND, '&Stand ID')
		selectMenu.AppendItem(sstd)
		selectMenu.AppendSeparator()
		sdig = wx.MenuItem(selectMenu, ID_SELECT_DIGITIZER, '&Digitizer Number')
		selectMenu.AppendItem(sdig)
		
		# Creating the menubar.
		menubar.Append(fileMenu, '&File')
		menubar.Append(colorMenu, '&Color Coding')
		menubar.Append(detailMenu, '&Details')
		menubar.Append(powerMenu, '&Average Power')
		menubar.Append(selectMenu, 'F&ind')
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
		self.Bind(wx.EVT_MENU, self.onSSMIF, id=ID_SSMIF)
		self.Bind(wx.EVT_MENU, self.onQuit, id=ID_QUIT)
		
		# Color menu events
		self.Bind(wx.EVT_MENU, self.onColor0, id=ID_COLOR_0)
		self.Bind(wx.EVT_MENU, self.onColor1, id=ID_COLOR_1)
		self.Bind(wx.EVT_MENU, self.onColor2, id=ID_COLOR_2)
		self.Bind(wx.EVT_MENU, self.onColor3, id=ID_COLOR_3)
		self.Bind(wx.EVT_MENU, self.onColor4, id=ID_COLOR_4)
		self.Bind(wx.EVT_MENU, self.onAdjust, id=ID_COLOR_ADJUST)
		
		# Detail menu events
		self.Bind(wx.EVT_MENU, self.onAntenna, id=ID_DETAIL_ANT)
		self.Bind(wx.EVT_MENU, self.onStand, id=ID_DETAIL_STAND)
		self.Bind(wx.EVT_MENU, self.onFEE, id=ID_DETAIL_FEE)
		self.Bind(wx.EVT_MENU, self.onCable, id=ID_DETAIL_CABLE)
		self.Bind(wx.EVT_MENU, self.onRFI, id=ID_DETAIL_RFI)
		self.Bind(wx.EVT_MENU, self.onStatus, id=ID_DETAIL_CHANGE_STATUS)
		
		# Power menu events
		self.Bind(wx.EVT_MENU, self.onAvgPower, id=ID_AVG_POWER)
		self.Bind(wx.EVT_MENU, self.onDataRange, id=ID_AVG_RANGE)
		self.Bind(wx.EVT_MENU, self.onAvgPowerSummary, id=ID_AVG_SUMMARY)
		
		# Select menu events
		self.Bind(wx.EVT_MENU, self.onSelectAntenna, id=ID_SELECT_ANTENNA)
		self.Bind(wx.EVT_MENU, self.onSelectStand, id=ID_SELECT_STAND)
		self.Bind(wx.EVT_MENU, self.onSelectDigitizer, id=ID_SELECT_DIGITIZER)
		
		# Make the images resizable
		self.Bind(wx.EVT_PAINT, self.resizePlots)
	
	def onOpen(self, event):
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
		
	def onSSMIF(self, event):
		"""
		Display the SSMIF antenna and FEE status codes.
		"""
		
		DisplaySSMIF(self)
		
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
		Set the stand field color coding to the antenna status.  Re-draw if 
		necessary.
		"""
		
		if self.data.color != 3:
			self.data.color = 3
			self.data.draw()
			if self.cAdjust is not None:
				try:
					self.cAdjust.Close()
				except:
					pass
			
	def onColor4(self, event):
		"""
		Set the stand field color coding to the estimates resonance point
		frequency in MHz.  Re-draw if necessary.
		"""
		
		if self.data.color != 4:
			self.data.color = 4
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
		
		if self.data.color != 3:
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
			
			fenDistA = numpy.zeros(4)
			k = 0
			for p1,p2 in zip([(-44.315,-72.522), (72.150,-44.277), (44.077,72.191), (-72.543,43.972)], [(72.150,-44.277), (44.077,72.191), (-72.543,43.972), (-44.315,-72.522)]):
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
			
			# Catch things outside the fence
			if abs(std.x) > 60 or abs(std.y) > 60:
				k = 0
				for p1 in [(-44.315,-72.522), (72.150,-44.277), (44.077,72.191), (-72.543,43.972)]:
					x1 = p1[0]
					y1 = p1[1]
					
					x3 = std.x
					y3 = std.y
					
					fenDistA[k] = numpy.sqrt( (x3-x1)**2 + (y3-y1)**2 )
					k += 1
				
			fenDist = fenDistA.min()
			
			outString = """Stand: %i

Coordinates:
x = %.3f m
y = %.3f m
z = %.3f m

Fence Distance: %.3f m
""" % (std.id, std.x, std.y, std.z, fenDist)
			
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
		
	def onStatus(self, event):
		"""
		Display the change antenna/FEE status dialog.
		"""
		
		if self.data.bestX > 0:
			StatusChangeDialog(self)
		
	def onAvgPower(self, event):
		"""
		Display the average power plots.
		"""
		
		if self.data.avgPower is not None and self.data.bestX > 0:
			AvgPowerDisplay(self)
			
	def onDataRange(self, event):
		"""
		Display the data range plots.
		"""
		
		if self.data.dataRange is not None and self.data.bestX > 0:
			DataRangeDisplay(self)
		
	def onAvgPowerSummary(self, event):
		"""
		Display a message box with the average power summary.
		"""
		
		if self.data is None:
			pass
		if self.data.avgPower is None:
			pass
		
		if self.data.bestX != -1:
			ant1 = self.data.antennas[self.data.bestX-1]
			dat1 = self.data.avgPower[self.data.bestX-1,:]
			ant2 = self.data.antennas[self.data.bestY-1]
			dat2 = self.data.avgPower[self.data.bestY-1,:]
		
			outString = """Antenna: %i
Polarization: %i

Global Mean:
%.2f +/- %.2f
Global Range:
%.2f to %.2f

Antenna: %i
Polarization: %i

Global Mean:
%.2f +/- %.2f
Global Range:
%.2f to %.2f
""" % (ant1.id, ant1.pol, dat1.mean(), dat1.std(), dat1.min(), dat1.max(), 
		ant2.id, ant2.pol, dat2.mean(), dat2.std(), dat2.min(), dat2.max())
		
			box = wx.MessageDialog(self, outString, "Average Power Details")
			box.ShowModal()
		else:
			pass
		
	def onSelectAntenna(self, event):
		"""
		Bring up a dialog box to find an antenna based on its ID number.
		"""
		
		box = SelectBox(self, mode='antenna')
		if box.ShowModal() == wx.ID_OK:
			antID = int(box.input.GetValue())
			if antID < 1 or antID > 520:
				pass
			elif self.data.antennas is None:
				pass
			else:
				for ant in self.data.antennas:
					if ant.id == antID:
						self.data.drawSpectrum(ant.stand.x, ant.stand.y)
						self.data.makeMark(ant.stand.x, ant.stand.y)
						break
				
		box.Destroy()
	
	def onSelectStand(self, event):
		"""
		Bring up a dialog box to find a stand based on its ID number.
		"""
		
		box = SelectBox(self, mode='stand')
		if box.ShowModal() == wx.ID_OK:
			stdID = int(box.input.GetValue())
			if stdID < 1 or stdID > 260:
				pass
			elif self.data.antennas is None:
				pass
			else:
				for ant in self.data.antennas:
					if ant.stand.id == stdID:
						self.data.drawSpectrum(ant.stand.x, ant.stand.y, preferStand=stdID)
						self.data.makeMark(ant.stand.x, ant.stand.y)
						break
				
		box.Destroy()
	
	def onSelectDigitizer(self, event):
		"""
		Bring up a dialog box to find a antenna associated with a particular 
		digitizer number.
		"""
		
		box = SelectBox(self, mode='digitizer')
		if box.ShowModal() == wx.ID_OK:
			digID = int(box.input.GetValue())
			if digID < 1 or digID > 520:
				pass
			elif self.data.antennas is None:
				pass
			else:
				for ant in self.data.antennas:
					if ant.digitizer == digID:
						self.data.drawSpectrum(ant.stand.x, ant.stand.y)
						self.data.makeMark(ant.stand.x, ant.stand.y)
						break
				
		box.Destroy()

	def resizePlots(self, event):
		w, h = self.GetSize()
		dpi = self.figure1.get_dpi()
		newW = 1.0*w/dpi
		newH1 = 1.0*(h/2-100)/dpi
		newH2 = 1.0*(h/2-75)/dpi
		self.figure1.set_size_inches((newW, newH1))
		self.figure1.canvas.draw()
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
		elif self.parent.data.color == 3:
			mode = 'Antenna Status'
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


class AvgPowerDisplay(wx.Frame):
	"""
	Window for displaying the average power with time for the selected stand.
	"""
	
	def __init__(self, parent):
		wx.Frame.__init__(self, parent, title='Time-Averaged Power', size=(400, 375))
		
		self.parent = parent
		
		self.initUI()
		self.initEvents()
		self.Show()
		
		self.initPlot()
		
	def __nextTen(self, value):
		"""
		Round a positive value to the next highest multiple of ten.
		"""
		
		return 10*numpy.ceil(value/10.0)
		
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
		Set all of the various events in the average power window.
		"""
		
		# Make the images resizable
		self.Bind(wx.EVT_PAINT, self.resizePlots)
		
	def initPlot(self):
		"""
		Populate the figure/canvas areas with a plot.  We only need to do this
		once for this type of window.
		"""
		
		avgPower = self.parent.data.avgPower
		bestX = self.parent.data.bestX
		bestY = self.parent.data.bestY
		
		if avgPower is None:
			return False
		if bestX < 1:
			return False
		
		self.figure.clf()
		self.ax1 = self.figure.gca()
		
		ant1 = self.parent.data.antennas[bestX-1]
		ant2 = self.parent.data.antennas[bestY-1]
		
		# Average power plot
		tScale = float(round(avgPower.shape[1] / 61.2244898))
		t = numpy.arange(0,avgPower.shape[1])/tScale + 0.5/tScale
		#self.ax1.plot(t, avgPower[bestX-1,:], label='Pol. %i' % ant1.pol)
		#self.ax1.plot(t, avgPower[bestY-1,:], label='Pol. %i' % ant2.pol)
		self.ax1.errorbar(t, avgPower[bestX-1,:], xerr=0.5/tScale, linestyle=' ', marker='+', label='Pol. %i' % ant1.pol, capsize=0)
		self.ax1.errorbar(t, avgPower[bestY-1,:], xerr=0.5/tScale, linestyle=' ', marker='+', label='Pol. %i' % ant2.pol, capsize=0)
		
		# Set ranges
		self.ax1.set_xlim([0, 61])
		self.ax1.set_ylim([0, self.__nextTen(avgPower.max())])
		
		# Labels
		self.ax1.set_title('Stand #%i' % ant1.stand.id)
		self.ax1.set_xlabel('Time [ms]')
		self.ax1.set_ylabel('Average Power [counts]')
		
		# Legend
		self.ax1.legend(loc=0)
		self.tScale = tScale
		
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
			
			try:
				tScale = self.tScale
				ap1 = self.parent.data.avgPower[self.parent.data.bestX-1,int(clickX*tScale)]
				ap2 = self.parent.data.avgPower[self.parent.data.bestY-1,int(clickX*tScale)]
				self.statusbar.SetStatusText("t=%.2f ms, X pol. Power=%.2f counts, Y pol. Power=%.2f" % (clickX, ap1, ap2))
			except IndexError:
				self.statusbar.SetStatusText("")
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


class DataRangeDisplay(wx.Frame):
	"""
	Window for displaying the time series mean, min, and maximum raw data 
	values.
	"""
	
	def __init__(self, parent):
		wx.Frame.__init__(self, parent, title='Range of Raw Data', size=(400, 375))
		
		self.parent = parent
		
		self.initUI()
		self.initEvents()
		self.Show()
		
		self.initPlot()
		
	def __nextTen(self, value):
		"""
		Round a positive value to the next highest multiple of ten.
		"""
		
		return 10*numpy.ceil(value/10.0)
		
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
		
		dataRange = self.parent.data.dataRange
		bestX = self.parent.data.bestX
		bestY = self.parent.data.bestY
		
		if dataRange is None:
			return False
		if bestX < 1:
			return False
		
		self.figure.clf()
		self.ax1 = self.figure.gca()
		
		ant1 = self.parent.data.antennas[bestX-1]
		ant2 = self.parent.data.antennas[bestY-1]
		
		# Data Range
		tScale = float(round(dataRange.shape[1] / 61.2244898))
		t = numpy.arange(0,dataRange.shape[1])/tScale + 0.5/tScale
		eb1 = numpy.zeros((2,t.size))
		eb1[0,:] = dataRange[bestX-1,:,1] - dataRange[bestX-1,:,0]
		eb1[1,:] = dataRange[bestX-1,:,2] - dataRange[bestX-1,:,1]
		self.ax1.errorbar(t, dataRange[bestX-1,:,1], xerr=0.5/tScale, yerr=eb1, linestyle=' ', marker='+', label='Pol. %i' % ant1.pol)
		eb2 = numpy.zeros((2,t.size))
		eb2[0,:] = dataRange[bestY-1,:,1] - dataRange[bestY-1,:,0]
		eb2[1,:] = dataRange[bestY-1,:,2] - dataRange[bestY-1,:,1]
		self.ax1.errorbar(t, dataRange[bestY-1,:,1], xerr=0.5/tScale, yerr=eb2, linestyle=' ', marker='+', label='Pol. %i' % ant2.pol)
		
		# Set ranges
		self.ax1.set_xlim([0, 61])
		self.ax1.set_ylim([-2048, 2047])
		
		# Labels
		self.ax1.set_title('Stand #%i' % ant1.stand.id)
		self.ax1.set_xlabel('Time [ms]')
		self.ax1.set_ylabel('Data Range [counts]')
		
		# Legend
		self.ax1.legend(loc=0)
		self.tScale = tScale
		
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
			
			try:
				tScale = self.tScale
				dr1 = self.parent.data.dataRange[self.parent.data.bestX-1,int(clickX*tScale),:]
				dr2 = self.parent.data.dataRange[self.parent.data.bestY-1,int(clickX*tScale),:]
				self.statusbar.SetStatusText("t=%.2f ms, X pol. Range: %+i to %+i, Y Pol. Range: %+i to %+i counts" % (clickX, dr1[0], dr1[2], dr2[0], dr2[2]))
			except IndexError:
				self.statusbar.SetStatusText("")
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


class SelectBox(wx.Dialog):
	"""
	Window for displaying the a simple dialog to find an antenna/stand/digitizer.
	"""
	
	def __init__(self, parent, mode='antenna'):
		wx.Dialog.__init__(self, parent, title='Find %s by ID' % mode.capitalize(), size=(200, 125))
		
		self.parent = parent
		self.mode = mode
		
		self.initUI()
		
	def initUI(self):
		"""
		Start the user interface.
		"""
		
		panel = wx.Panel(self, -1)
		vbox = wx.BoxSizer(wx.VERTICAL)
		
		wx.StaticBox(panel, -1, '%s ID' % self.mode.capitalize(), (5, 5), (190, 75))
		self.input = wx.TextCtrl(panel, -1, '', (15, 30))
		if self.mode == 'stand':
			wx.StaticText(panel, -1, 'Limits: 1 - 260, inclusive', (15, 60))
		else:
			wx.StaticText(panel, -1, 'Limits: 1 - 520, inclusive', (15, 60))

		hbox = wx.BoxSizer(wx.HORIZONTAL)
		okButton = wx.Button(self, wx.ID_OK, 'Ok', size=(70, 30))
		closeButton = wx.Button(self, wx.ID_CANCEL, 'Close', size=(70, 30))
		hbox.Add(okButton, 1)
		hbox.Add(closeButton, 1, wx.LEFT, 5)

		vbox.Add(panel)
		vbox.Add(hbox, 1, wx.ALIGN_CENTER | wx.TOP | wx.BOTTOM, 10)

		self.SetSizer(vbox)


STATUS_CHANGE_OK = 201
STATUS_CHANGE_CANCEL = 202

class StatusChangeDialog(wx.Frame):
	"""
	Window for chaning the status of an antenna or FEE.
	"""
	
	def __init__ (self, parent):	
		wx.Frame.__init__(self, parent, title='Change Antenna/FEE Status', size=(450, 225))
		
		self.parent = parent
		
		self.initUI()
		self.initEvents()
		self.Show()
		
	def initUI(self):
		bestX = self.parent.data.bestX
		bestY = self.parent.data.bestY
		
		ant1 = self.parent.data.antennas[bestX-1]
		ant2 = self.parent.data.antennas[bestY-1]
		
		row = 0
		panel = wx.Panel(self)
		sizer = wx.GridBagSizer(5, 5)
		
		std = wx.StaticText(panel, label='Stand #%i' % ant1.stand.id)
		sizer.Add(std, pos=(row+0, 0), span=(1,4), flag=wx.EXPAND|wx.LEFT|wx.RIGHT, border=5)
		row += 1
		
		polX = wx.StaticText(panel, label='X Pol.')
		sizer.Add(polX, pos=(row+0, 0), span=(1,4), flag=wx.EXPAND|wx.LEFT|wx.RIGHT, border=5)
		row += 1
		
		asX = wx.StaticText(panel, label='Antenna Status')
		sizer.Add(asX, pos=(row+0, 0), span=(1,1), flag=wx.EXPAND|wx.LEFT|wx.RIGHT, border=5)
		self.as10 = wx.RadioButton(panel, -1, 'Not Installed', style=wx.RB_GROUP)
		self.as11 = wx.RadioButton(panel, -1, 'Bad')
		self.as12 = wx.RadioButton(panel, -1, 'Suspect')
		self.as13 = wx.RadioButton(panel, -1, 'Good')
		sizer.Add(self.as10, pos=(row+0, 1), span=(1,1))
		sizer.Add(self.as11, pos=(row+0, 2), span=(1,1))
		sizer.Add(self.as12, pos=(row+0, 3), span=(1,1))
		sizer.Add(self.as13, pos=(row+0, 4), span=(1,1))
		row += 1
		
		polY = wx.StaticText(panel, label='Y Pol.')
		sizer.Add(polY, pos=(row+0, 0), span=(1,4), flag=wx.EXPAND|wx.LEFT|wx.RIGHT, border=5)
		row += 1
		
		asY = wx.StaticText(panel, label='Antenna Status')
		sizer.Add(asY, pos=(row+0, 0), span=(1,1), flag=wx.EXPAND|wx.LEFT|wx.RIGHT, border=5)
		self.as20 = wx.RadioButton(panel, -1, 'Not Installed', style=wx.RB_GROUP)
		self.as21 = wx.RadioButton(panel, -1, 'Bad')
		self.as22 = wx.RadioButton(panel, -1, 'Suspect')
		self.as23 = wx.RadioButton(panel, -1, 'Good')
		sizer.Add(self.as20, pos=(row+0, 1), span=(1,1))
		sizer.Add(self.as21, pos=(row+0, 2), span=(1,1))
		sizer.Add(self.as22, pos=(row+0, 3), span=(1,1))
		sizer.Add(self.as23, pos=(row+0, 4), span=(1,1))
		row += 1
		
		line = wx.StaticLine(panel)
		sizer.Add(line, pos=(row+0, 0), span=(1, 5), flag=wx.EXPAND|wx.BOTTOM, border=10)
		row += 1
		
		fs =  wx.StaticText(panel, label='FEE Status')
		sizer.Add(fs, pos=(row+0, 0), span=(1,1), flag=wx.EXPAND|wx.LEFT|wx.RIGHT, border=5)
		self.fs0 = wx.RadioButton(panel, -1, 'Not Installed', style=wx.RB_GROUP)
		self.fs1 = wx.RadioButton(panel, -1, 'Bad')
		self.fs2 = wx.RadioButton(panel, -1, 'Suspect')
		self.fs3 = wx.RadioButton(panel, -1, 'Good')
		sizer.Add(self.fs0, pos=(row+0, 1), span=(1,1))
		sizer.Add(self.fs1, pos=(row+0, 2), span=(1,1))
		sizer.Add(self.fs2, pos=(row+0, 3), span=(1,1))
		sizer.Add(self.fs3, pos=(row+0, 4), span=(1,1))
		row += 1
		
		line = wx.StaticLine(panel)
		sizer.Add(line, pos=(row+0, 0), span=(1, 5), flag=wx.EXPAND|wx.BOTTOM, border=10)
		row += 1
		
		#
		# Buttons
		#
		ok = wx.Button(panel, STATUS_CHANGE_OK, 'Ok', size=(56, 28))
		cancel = wx.Button(panel, STATUS_CHANGE_CANCEL, 'Cancel', size=(56, 28))
		sizer.Add(ok, pos=(row+0, 3), flag=wx.RIGHT|wx.BOTTOM, border=5)
		sizer.Add(cancel, pos=(row+0, 4), flag=wx.RIGHT|wx.BOTTOM, border=5)
		
		panel.SetSizerAndFit(sizer)
		
		#
		# Set current values
		#
		
		## Antenna 1
		if ant1.status == 0:
			self.as10.SetValue(True)
			self.as11.SetValue(False)
			self.as12.SetValue(False)
			self.as13.SetValue(False)
		elif ant1.status == 1:
			self.as10.SetValue(False)
			self.as11.SetValue(True)
			self.as12.SetValue(False)
			self.as13.SetValue(False)
		elif ant1.status == 2:
			self.as10.SetValue(False)
			self.as11.SetValue(False)
			self.as12.SetValue(True)
			self.as13.SetValue(False)
		else:
			self.as10.SetValue(False)
			self.as11.SetValue(False)
			self.as12.SetValue(False)
			self.as13.SetValue(True)
		
		## Antenna 2
		if ant2.status == 0:
			self.as20.SetValue(True)
			self.as21.SetValue(False)
			self.as22.SetValue(False)
			self.as23.SetValue(False)
		elif ant2.status == 1:
			self.as20.SetValue(False)
			self.as21.SetValue(True)
			self.as22.SetValue(False)
			self.as23.SetValue(False)
		elif ant2.status == 2:
			self.as20.SetValue(False)
			self.as21.SetValue(False)
			self.as22.SetValue(True)
			self.as23.SetValue(False)
		else:
			self.as20.SetValue(False)
			self.as21.SetValue(False)
			self.as22.SetValue(False)
			self.as23.SetValue(True)
		
		## FEE
		if ant1.fee.status == 0:
			self.fs0.SetValue(True)
			self.fs1.SetValue(False)
			self.fs2.SetValue(False)
			self.fs3.SetValue(False)
		elif ant1.fee.status == 1:
			self.fs0.SetValue(False)
			self.fs1.SetValue(True)
			self.fs2.SetValue(False)
			self.fs3.SetValue(False)
		elif ant1.fee.status == 2:
			self.fs0.SetValue(False)
			self.fs1.SetValue(False)
			self.fs2.SetValue(True)
			self.fs3.SetValue(False)
		else:
			self.fs0.SetValue(False)
			self.fs1.SetValue(False)
			self.fs2.SetValue(False)
			self.fs3.SetValue(True)
		
	def initEvents(self):
		self.Bind(wx.EVT_BUTTON, self.onOk, id=STATUS_CHANGE_OK)
		self.Bind(wx.EVT_BUTTON, self.onCancel, id=STATUS_CHANGE_CANCEL)
		
	def onOk(self, event):
		# Antenna 1
		if self.as10.GetValue():
			self.parent.data.antennas[self.parent.data.bestX-1].status = 0
		elif self.as11.GetValue():
			self.parent.data.antennas[self.parent.data.bestX-1].status = 1
		elif self.as12.GetValue():
			self.parent.data.antennas[self.parent.data.bestX-1].status = 2
		else:
			self.parent.data.antennas[self.parent.data.bestX-1].status = 3
			
		# Antenna 2
		if self.as20.GetValue():
			self.parent.data.antennas[self.parent.data.bestY-1].status = 0
		elif self.as21.GetValue():
			self.parent.data.antennas[self.parent.data.bestY-1].status = 1
		elif self.as22.GetValue():
			self.parent.data.antennas[self.parent.data.bestY-1].status = 2
		else:
			self.parent.data.antennas[self.parent.data.bestY-1].status = 3
		
		# FEE
		if self.fs0.GetValue():
			self.parent.data.antennas[self.parent.data.bestX-1].fee.status = 0
			self.parent.data.antennas[self.parent.data.bestY-1].fee.status = 0
		elif self.fs1.GetValue():
			self.parent.data.antennas[self.parent.data.bestX-1].fee.status = 1
			self.parent.data.antennas[self.parent.data.bestY-1].fee.status = 1
		elif self.fs2.GetValue():
			self.parent.data.antennas[self.parent.data.bestX-1].fee.status = 2
			self.parent.data.antennas[self.parent.data.bestY-1].fee.status = 2
		else:
			self.parent.data.antennas[self.parent.data.bestX-1].fee.status = 3
			self.parent.data.antennas[self.parent.data.bestY-1].fee.status = 3
			
		# Refresh if we are in the antenna status color coding
		if self.parent.data.color == 3:
			self.parent.data.draw()
		
		self.Close()
		
	def onCancel(self, event):
		self.Close()


SSMIF_OK = 301

class DisplaySSMIF(wx.Frame):
	"""
	Text display window for printing out the new SSMIF entries for antenna status.
	FEE status is currently not supported because of a limitation in the LSL SSMIF
	parser.
	"""
	
	def __init__ (self, parent):	
		wx.Frame.__init__(self, parent, title='SSMIF Status Codes', size=(600, 600))
		
		self.parent = parent
		
		self.initUI()
		self.initEvents()
		self.generateText()
		self.Show()
		
	def initUI(self):
		vbox = wx.BoxSizer(wx.VERTICAL)
		self.textCtrl = wx.TextCtrl(self, -1, "", style=wx.TE_MULTILINE|wx.TE_READONLY, size=(600,500))
		vbox.Add(self.textCtrl, 1, wx.LEFT|wx.RIGHT|wx.BOTTOM|wx.TOP|wx.EXPAND, border=5)
		
		
		ok = wx.Button(self, SSMIF_OK, 'Ok', size=(56, 28))
		vbox.Add(ok, 0, flag=wx.ALIGN_RIGHT|wx.LEFT|wx.RIGHT|wx.BOTTOM, border=5)
		
		self.SetSizer(vbox)
		self.SetAutoLayout(1)
		vbox.Fit(self)
			
	def initEvents(self):
		self.Bind(wx.EVT_BUTTON, self.onOk,  id=SSMIF_OK)
		
	def generateText(self):
		def sortAnts(x, y):
			"""
			Small function to re-sort a list of Antenna instances by antenna number.
			"""
			
			if x.id > y.id:
				return 1
			elif x.id < y.id:
				return -1
			else:
				return 0
		# Sort antennas by antenna number
		ants = sorted(self.parent.data.antennas, cmp=sortAnts)
		
		#
		# Antenna status codes
		#
		
		self.textCtrl.AppendText('# -----------------------------\n# --- Antenna Status ---\n# -----------------------------\n# Status codes 0-3 summarized defined at end of this document (and in MCS0031)\n# This refers to the *antenna*, not the FEE or some combination of the two.\n# This will be set to 3 ("OK") for any antenna n <= 2*N_STD not identified.\n# *** ANT_STAT[antenna_id] goes here:\n')
		for ant in ants:
			if ant.status == 3:
				continue
			
			if ant.id < 10:
				self.textCtrl.AppendText('ANT_STAT[%i]    %i\n' % (ant.id, ant.status))
			elif ant.id < 100:
				self.textCtrl.AppendText('ANT_STAT[%i]   %i\n' % (ant.id, ant.status))
			else:
				self.textCtrl.AppendText('ANT_STAT[%i]  %i\n' % (ant.id, ant.status))
		self.textCtrl.AppendText('\n\n')
		
		#
		# FEE status codes
		#
		
		#self.textCtrl.AppendText('# ----------------------\n# --- FEE Status ---\n# ----------------------\n# Status codes 0-3 summarized defined at end of this document (and in MCS0031)\n# This will be set to 3 ("OK") for any FEE #\'s <= N_FEE not identified\n# *** FEE_STAT[fee_id] goes here:\n')
		#for ant in ants:
			#if ant.pol == 0:
				#self.textCtrl.AppendText('FEE_STAT[%s]  %i\n' % (ant.fee.id, ant.fee.status))
		#self.textCtrl.AppendText('\n')
		
		self.textCtrl.ShowPosition(0)
		
	def onOk(self, event):
		self.Close()


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
