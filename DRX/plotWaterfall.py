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
		self.crange = [timeBlocks.min()*1.05, timeBlocks.max()*1.05]
		
		self.mask = numpy.zeros((timeBlocks.shape[1], freq.shape[0]))
		
		self.connect()
		wx.EndBusyCursor()
		self.frame.enableButtons()
	
	def draw(self):
		"""Draw the waterfall diagram"""
		
		wx.BeginBusyCursor()
		self.frame.disableButtons()
		self.frame.setTextRange()
		
		freq = self.data[0]
		spec = self.data[1][:,self.index,:]
		spec = numpy.ma.array(spec, mask=numpy.zeros_like(spec))
		for i in xrange(spec.shape[0]):
			spec.mask[i,:] = self.mask[self.index,:]
		
		self.frame.figure1a.clf()
		ax = self.frame.figure1a.gca()
		m = ax.imshow(spec, 
					extent=(freq[0]/1e6, freq[-1]/1e6, 0, self.intTime*(spec.shape[0]-1)), 
					vmin=self.crange[0], vmax=self.crange[1])
		cm = self.frame.figure1a.colorbar(m, ax=ax)
		cm.ax.set_ylabel('PSD [arb. dB]')
		ax.axis('auto')
		ax.set_xlabel('Frequency [MHz]')
		ax.set_ylabel('Time [s]')
		ax.set_title('Tuning %i, Pol. %s' % (self.index/2+1, 'Y' if self.index %2 else 'X'))
		
		try:
			if self.oldLine1a is not None:
				ax.lines.extend(self.oldLine1a)
		except AttributeError:
			pass
		
		self.frame.canvas1a.draw()
		
		tp = 10.0**(spec/10.0)
		tp = tp.sum(axis=1)
		tp = numpy.log10(tp)*10
		
		self.frame.figure1b.clf()
		ax = self.frame.figure1b.gca()
		ax.plot(tp, numpy.arange(spec.shape[0])*self.intTime)
		ax.set_ylim([0, self.intTime*(spec.shape[0]-1)])
		ax.set_xlabel('Total Power [arb. dB]')
		ax.set_ylabel('Time [s]')
		
		try:
			if self.oldLine1b is not None:
				ax.lines.extend(self.oldLine1b)
		except AttributeError:
			pass
		
		self.frame.canvas1b.draw()
		
		wx.EndBusyCursor()
		self.frame.enableButtons()
	
	def drawSpectrum(self, clickY):
		"""Get the spectrum at a particular point in time."""

		freq = self.data[0]
		spec = self.data[1][int(round(clickY)),self.index,:]
		spec = numpy.ma.array(spec, mask=self.mask[self.index,:])
		medianSpec = numpy.median(self.data[1][:,self.index,:], axis=0)
		
		self.frame.figure2.clf()
		ax2 = self.frame.figure2.gca()
		ax2.plot(freq/1e6, spec, label='Current')
		ax2.plot(freq/1e6, medianSpec, label='Median', alpha=0.5)
		ax2.legend(loc=0)
		ax2.set_ylim(self.crange)
		ax2.set_xlabel('Frequency [MHz]')
		ax2.set_ylabel('PSD [arb. dB]')
		
		self.frame.canvas2.draw()
		self.spectrumClick = clickY
	
	def makeMark(self, clickY):
		ax = self.frame.figure1a.gca()
		
		try:
			if self.oldLine1a is not None:
				del ax.lines[-1]
		except AttributeError:
			pass
		
		freq = self.data[0]
			
		oldXSize = ax.get_xlim()
		oldYSize = ax.get_ylim()
		
		self.oldLine1a = ax.plot(freq/1e6, freq*0+round(clickY), color='red')
		ax.set_xlim(oldXSize)
		ax.set_ylim(oldYSize)
		
		self.frame.canvas1a.draw()
		
		###
		
		ax = self.frame.figure1b.gca()
		
		try:
			if self.oldLine1b is not None:
				del ax.lines[-1]
		except AttributeError:
			pass
			
		oldXSize = ax.get_xlim()
		oldYSize = ax.get_ylim()
		
		self.oldLine1b = ax.plot(oldXSize, [round(clickY)]*2, color='red')
		ax.set_xlim(oldXSize)
		ax.set_ylim(oldYSize)
		
		self.frame.canvas1b.draw()
	
	def connect(self):
		'connect to all the events we need'
		
		self.cidpress1a = self.frame.figure1a.canvas.mpl_connect('button_press_event', self.on_press)
		self.cidpress1b = self.frame.figure1b.canvas.mpl_connect('button_press_event', self.on_press)
		self.cidpress2 = self.frame.figure2.canvas.mpl_connect('button_press_event', self.on_press2)
		self.cidmotion = self.frame.figure1a.canvas.mpl_connect('motion_notify_event', self.on_motion)
	
	def on_press(self, event):
		'on button press we will see if the mouse is over us and store some data'
		
		if event.inaxes:
			clickX = event.xdata
			clickY = event.ydata / self.intTime
			
			freq = self.data[0]
			
			self.draw()
			self.drawSpectrum(clickY)
			self.makeMark(clickY*self.intTime)
			
	def on_press2(self, event):
		if event.inaxes:
			clickX = event.xdata*1e6
			clickY = event.ydata
			
			freq = self.data[0]
			best = numpy.where( numpy.abs(freq-clickX) == numpy.abs(freq-clickX).min() )[0][0]
			
			if event.button == 3:
				self.mask[self.index, best] = 1
			elif event.button == 2:
				self.mask[self.index, best] = 0
			else:
				pass
			self.drawSpectrum(self.spectrumClick)
			
	def on_motion(self, event):
		if event.inaxes:
			clickX = event.xdata
			clickY = event.ydata
			clickYp = clickY / self.intTime
			
			dataX = numpy.where(numpy.abs(clickX-self.data[0]/1e6) == (numpy.abs(clickX-self.data[0]/1e6).min()))[0][0]
			
			value = self.data[1][int(round(clickYp)), self.index, int(round(dataX))]
			self.frame.statusbar.SetStatusText("f=%.4f MHz, t=%.4f s, p=%.2f dB" % (clickX, clickY, value))
		else:
			self.frame.statusbar.SetStatusText("")
			
	
	def disconnect(self):
		'disconnect all the stored connection ids'
		
		self.frame.figure1a.canvas.mpl_disconnect(self.cidpress1a)
		self.frame.figure1b.canvas.mpl_disconnect(self.cidpress1b)
		self.frame.figure2.canvas.mpl_disconnect(self.cidpress2)
		self.frame.figure1a.canvas.mpl_disconnect(self.cidmotion)


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
		hbox1 = wx.BoxSizer(wx.HORIZONTAL)
		self.figure1a = Figure()
		self.canvas1a = FigureCanvasWxAgg(panel1, -1, self.figure1a)
		self.figure1b = Figure()
		self.canvas1b = FigureCanvasWxAgg(panel1, -1, self.figure1b)
		
		hbox1.Add(self.canvas1a, 1, wx.EXPAND)
		hbox1.Add(self.canvas1b, 1, wx.EXPAND)
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
		#self.Bind(wx.EVT_PAINT, self.resizePlots)
	
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
		
	def setTextRange(self):
		low, hig = self.data.crange
		self.maxText.SetLabel("%+.2f" % hig)
		self.minText.SetLabel("%+.2f" % low)
		
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
