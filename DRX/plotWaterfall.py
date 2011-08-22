#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Given the NPZ output of drxWaterfall, plot it in an interative way.

$Rev$
$LastChangedBy$
$LastChangedDate$
"""

import os
import wx
import sys
import numpy

import matplotlib
matplotlib.use('WXAgg')
matplotlib.interactive(True)
from matplotlib.backends.backend_wxagg import NavigationToolbar2WxAgg, FigureCanvasWxAgg
from matplotlib.figure import Figure


class Waterfall(object):
	def __init__(self, frame, data=None):
		self.frame = frame
		self.data = data
		self.press = None
		
		self.index = 0
		self.ax = None
		self.ax2 = None
		
		self.spectrumClick = None
		
	def loadData(self, filename):
		"""Load in data from an NPZ file."""
		
		wx.BeginBusyCursor()
		self.frame.disableButtons()
		
		self.filename = filename
		dataDictionary = numpy.load(filename)
		
		freq = dataDictionary['freq']
		LFFT = freq.size
		timeBlocks = numpy.log10(dataDictionary['spec'])*10.0
		integrationTime = dataDictionary['tInt']
		
		self.orig = dataDictionary
		self.data = (freq, timeBlocks)
		self.intTime = integrationTime
		self.xrange = [freq.min(), freq.max()]
		self.yrange = [0, timeBlocks.shape[0]]
		self.crange = [-10.0, 20.0]
		
		self.connect()
		wx.EndBusyCursor()
		self.frame.enableButtons()
	
	def draw(self):
		"""Draw the waterfall diagram"""
		
		wx.BeginBusyCursor()
		self.frame.disableButtons()
		self.frame.setTextRange()
		
		self.frame.figure.clf()
		ax = self.frame.figure.gca()
		m = ax.imshow(self.data[1][:,self.index,:], 
					extent=(self.data[0][0]/1e6, self.data[0][-1]/1e6, 0, self.intTime*(self.data[1].shape[0]-1)), 
					vmin=self.crange[0], vmax=self.crange[1])
		cm = self.frame.figure.colorbar(m, ax=ax)
		cm.ax.set_ylabel('PSD [arb. dB]')
		ax.axis('auto')
		ax.set_xlabel('Frequency [MHz]')
		ax.set_ylabel('Time [s]')
		ax.set_title('Tuning %i, Pol. %s' % (self.index/2+1, 'Y' if self.index %2 else 'X'))
		
		try:
			if self.oldLine is not None:
				ax.lines.extend(self.oldLine)
		except AttributeError:
			pass
		
		self.frame.canvas.draw()
		
		wx.EndBusyCursor()
		self.frame.enableButtons()
	
	def drawSpectrum(self, clickY):
		"""Get the spectrum at a particular point in time."""
		
		freq = self.data[0]
		spec = self.data[1][int(round(clickY)),self.index,:]
		
		self.frame.figure2.clf()
		ax2 = self.frame.figure2.gca()
		ax2.plot(freq/1e6, spec)
		ax2.set_xlabel('Frequency [MHz]')
		ax2.set_ylabel('PSD [arb. dB]')
		
		self.frame.canvas2.draw()
		self.spectrumClick = clickY
	
	def makeMark(self, clickY):
		ax = self.frame.figure.gca()
		
		try:
			if self.oldLine is not None:
				del ax.lines[-1]
		except AttributeError:
			pass
		
		freq = self.data[0]
			
		oldXSize = ax.get_xlim()
		oldYSize = ax.get_ylim()
		
		self.oldLine = ax.plot(freq/1e6, freq*0+round(clickY), color='red')
		ax.set_xlim(oldXSize)
		ax.set_ylim(oldYSize)
		
		self.frame.canvas.draw()
	
	def connect(self):
		'connect to all the events we need'
		
		self.cidpress = self.frame.figure.canvas.mpl_connect('button_press_event', self.on_press)
		#self.cidrelease = self.figure.canvas.mpl_connect('button_release_event', self.on_release)
		self.cidmotion = self.frame.figure.canvas.mpl_connect('motion_notify_event', self.on_motion)
	
	def on_press(self, event):
		'on button press we will see if the mouse is over us and store some data'
		
		if event.inaxes and self.frame.toolbar.mode == '':
			clickX = event.xdata
			clickY = event.ydata
			
			freq = self.data[0]
			
			self.drawSpectrum(clickY)
			self.makeMark(clickY)
			
	def on_motion(self, event):
		if event.inaxes:
			clickX = event.xdata
			clickY = event.ydata
			
			dataX = numpy.where(numpy.abs(clickX-self.data[0]/1e6) == (numpy.abs(clickX-self.data[0]/1e6).min()))[0][0]
			
			value = self.data[1][int(round(clickY)), self.index, int(round(dataX))]
			self.frame.statusbar.SetStatusText("f=%.4f MHz, t=%.4f s, p=%.2f dB" % (clickX, clickY, value))
		else:
			self.frame.statusbar.SetStatusText("")
			
	
	def disconnect(self):
		'disconnect all the stored connection ids'
		
		self.frame.figure.canvas.mpl_disconnect(self.cidpress)
		#self.figure.canvas.mpl_disconnect(self.cidrelease)
		self.frame.figure.canvas.mpl_disconnect(self.cidmotion)


ID_OPEN = 10
ID_QUIT = 11
ID_TUNING1_X = 20
ID_TUNING1_Y = 21
ID_TUNING2_X = 22
ID_TUNING2_Y = 23

class MainWindow(wx.Frame):
	def __init__(self, parent, id):
		self.dirname = ''
		self.data = None
		
		wx.Frame.__init__(self, parent, id, title="DRX Waterfall Viewer", size=(600,800))
		
		self.initUI()
		self.initEvents()
		self.Show()
		
	def initUI(self):
		self.statusbar = self.CreateStatusBar() # A Statusbar in the bottom of the window
		
		font = wx.SystemSettings_GetFont(wx.SYS_SYSTEM_FONT)
		font.SetPointSize(10)
		
		menuBar = wx.MenuBar()
		
		fileMenu = wx.Menu()
		open = wx.MenuItem(fileMenu, ID_OPEN, "&Open")
		fileMenu.AppendItem(open)
		fileMenu.AppendSeparator()
		exit = wx.MenuItem(fileMenu, ID_QUIT, "E&xit")
		fileMenu.AppendItem(exit)
		
		dataMenu = wx.Menu()
		dataMenu.AppendRadioItem(ID_TUNING1_X, 'Tuning 1, Pol. X')
		dataMenu.AppendRadioItem(ID_TUNING1_Y, 'Tuning 1, Pol. Y')
		dataMenu.AppendSeparator()
		dataMenu.AppendRadioItem(ID_TUNING2_X, 'Tuning 2, Pol. X')
		dataMenu.AppendRadioItem(ID_TUNING2_Y, 'Tuning 2, Pol. Y')

		# Creating the menubar.
		menuBar.Append(fileMenu,"&File") # Adding the "filemenu" to the MenuBar
		menuBar.Append(dataMenu, "&Data")
		self.SetMenuBar(menuBar)  # Adding the MenuBar to the Frame content.
		
		vbox = wx.BoxSizer(wx.VERTICAL)
		
		# Add waterfall plot
		panel1 = wx.Panel(self, -1)
		hbox1 = wx.BoxSizer(wx.VERTICAL)
		self.figure = Figure()
		self.canvas = FigureCanvasWxAgg(panel1, -1, self.figure)
		self.toolbar = NavigationToolbar2WxAgg(self.canvas)
		self.toolbar.Realize()
		hbox1.Add(self.toolbar, 0, wx.LEFT | wx.FIXED_MINSIZE)
		hbox1.Add(self.canvas, 1, wx.EXPAND)
		panel1.SetSizer(hbox1)
		vbox.Add(panel1, 1, wx.EXPAND)
		
		# Add a contrast control panel
		self.plotButtons = []
		panel2 = wx.Panel(self, -1)
		hbox2 = wx.BoxSizer(wx.HORIZONTAL)
		## Range
		sizerRange = wx.StaticBoxSizer(wx.StaticBox(panel2, -1, 'Range'), orient=wx.VERTICAL)
		self.plotButtons.append(wx.Button(panel2, -1, "Decrease", size=(120, -1)))
		self.Bind(wx.EVT_BUTTON, self.decRange, self.plotButtons[-1])
		sizerRange.Add(self.plotButtons[-1])
		self.plotButtons.append(wx.Button(panel2, -1, "Increase", size=(120, -1)))
		self.Bind(wx.EVT_BUTTON, self.incRange, self.plotButtons[-1])
		sizerRange.Add(self.plotButtons[-1])
		hbox2.Add(sizerRange, 0, wx.EXPAND)
		## Clip
		sizerClip = wx.StaticBoxSizer(wx.StaticBox(panel2, -1, 'Clip'), orient=wx.VERTICAL)
		self.plotButtons.append(wx.Button(panel2, -1, "Decrease", size=(120,-1)))
		self.Bind(wx.EVT_BUTTON, self.decClip, self.plotButtons[-1])
		sizerClip.Add(self.plotButtons[-1])
		self.plotButtons.append(wx.Button(panel2, -1, "Increase", size=(120,-1)))
		self.Bind(wx.EVT_BUTTON, self.incClip, self.plotButtons[-1])
		sizerClip.Add(self.plotButtons[-1])
		hbox2.Add(sizerClip, 0, wx.EXPAND)
		## Display range
		sizerMap = wx.StaticBoxSizer(wx.StaticBox(panel2, -1, 'Color Map'), orient=wx.VERTICAL)
		gridMap = wx.GridSizer(2, 3)
		gridMap.Add(wx.StaticText(panel2, id=-1, label='Max (red):  '))
		self.maxText = wx.StaticText(panel2, id=-1, label='--', style=wx.ALIGN_RIGHT)
		gridMap.Add(self.maxText)
		gridMap.Add(wx.StaticText(panel2, id=-1, label='dB'))
		gridMap.Add(wx.StaticText(panel2, id=-1, label='Min (blue):  '))
		self.minText = wx.StaticText(panel2, id=-1, label='--', style=wx.ALIGN_RIGHT)
		gridMap.Add(self.minText)
		gridMap.Add(wx.StaticText(panel2, id=-1, label='dB'))
		sizerMap.Add(gridMap)
		hbox2.Add(sizerMap, 0, wx.EXPAND)
		
		panel2.SetSizer(hbox2)
		vbox.Add(panel2, 0, wx.LEFT)
		
		# Add a spectrum plot
		panel3 = wx.Panel(self, -1)
		hbox3 = wx.BoxSizer(wx.HORIZONTAL)
		self.figure2 = Figure()
		self.canvas2 = FigureCanvasWxAgg(panel3, -1, self.figure2)
		hbox3.Add(self.canvas2, 1)
		panel3.SetSizer(hbox3)
		vbox.Add(panel3, 1, wx.EXPAND)
		
		# Use some sizers to see layout options
		self.SetSizer(vbox)
		self.SetAutoLayout(1)
		vbox.Fit(self)
		
		#Layout sizers
		self.Show()
		
	def initEvents(self):
		self.Bind(wx.EVT_MENU, self.onOpen, id=ID_OPEN)
		self.Bind(wx.EVT_MENU, self.onExit, id=ID_QUIT)
		
		self.Bind(wx.EVT_MENU, self.onTuning1X, id=ID_TUNING1_X)
		self.Bind(wx.EVT_MENU, self.onTuning1Y, id=ID_TUNING1_Y)
		self.Bind(wx.EVT_MENU, self.onTuning2X, id=ID_TUNING2_X)
		self.Bind(wx.EVT_MENU, self.onTuning2Y, id=ID_TUNING2_Y)
		
		# Make the images resizable
		self.Bind(wx.EVT_PAINT, self.resizePlots)
	
	def onOpen(self,e):
		""" Open a file"""
		dlg = wx.FileDialog(self, "Choose a file", self.dirname, "", "*.*", wx.OPEN)
		if dlg.ShowModal() == wx.ID_OK:
			self.filename = dlg.GetFilename()
			self.dirname = dlg.GetDirectory()
			self.data = Waterfall(self)
			self.data.loadData(os.path.join(self.dirname, self.filename))
			self.data.draw()
		dlg.Destroy()
		
	def onExit(self, event):
		self.Close(True)
		
	def onTuning1X(self, event):
		self.data.index = 0
		self.data.draw()
		if self.data.spectrumClick is not None:
			self.data.drawSpectrum(self.data.spectrumClick)
		
	def onTuning1Y(self, event):
		self.data.index = 1
		self.data.draw()
		if self.data.spectrumClick is not None:
			self.data.drawSpectrum(self.data.spectrumClick)
		
	def onTuning2X(self, event):
		self.data.index = 2
		self.data.draw()
		if self.data.spectrumClick is not None:
			self.data.drawSpectrum(self.data.spectrumClick)
		
	def onTuning2Y(self, event):
		self.data.index = 3
		self.data.draw()
		if self.data.spectrumClick is not None:
			self.data.drawSpectrum(self.data.spectrumClick)

	def decRange(self, event):
		if self.data is not None:
			low, hig = self.data.crange
			low += 1
			hig -=1
			self.data.crange = [low, hig]
			self.data.draw()
		
	def incRange(self, event):
		if self.data is not None:
			low, hig = self.data.crange
			low -= 1
			hig +=1
			self.data.crange = [low, hig]
			self.data.draw()
			
	def decClip(self, event):
		if self.data is not None:
			low, hig = self.data.crange
			low -= 1
			hig -= 1
			self.data.crange = [low, hig]
			self.data.draw()
		
	def incClip(self, event):
		if self.data is not None:
			low, hig = self.data.crange
			low += 1
			hig += 1
			self.data.crange = [low, hig]
			self.data.draw()
			
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
		
	def setTextRange(self):
		low, hig = self.data.crange
		self.maxText.SetLabel("%+.2f" % hig)
		self.minText.SetLabel("%+.2f" % low)

	def GetToolBar(self):
		# You will need to override GetToolBar if you are using an 
		# unmanaged toolbar in your frame
		return self.toolbar
		
	def disableButtons(self):
		for button in self.plotButtons:
			button.Enable(False)
			
	def enableButtons(self):
		for button in self.plotButtons:
			button.Enable(True)

def main(args):
	app = wx.App(0)
	frame = MainWindow(None, -1)
	if len(args) == 1:
		frame.data = Waterfall(frame)
		frame.data.loadData(args[0])
		frame.data.draw()
	app.MainLoop()


if __name__ == "__main__":
	main(sys.argv[1:])
