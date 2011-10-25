#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import numpy
from lsl.common.stations import lwa1
from lsl.statistics import robust
from matplotlib import pyplot as plt


def clean180Jumps(x):
	factor = 0.70
	
	out = 1.0*x
	for i in xrange(1,out.size):
		if out[i-1] - out[i] < -factor*2*numpy.pi:
			out[i] -= 2*numpy.pi
		elif out[i-1] - out[i] > factor*2*numpy.pi:
			out[i] += 2*numpy.pi
			
	return out


class FitRegion:
	def __init__(self, ax1, ax2, freq, phase, extra=None):
		self.freq  = freq
		self.phase = phase
		self.corr  = None
		self.extra = extra
		
		self.nClick1 = 0
		self.x11 = None
		self.x12 = None
		self.ax1 = ax1
		self.ax2 = ax2
		self.cid1 = self.ax1.figure.canvas.mpl_connect('button_press_event', self)

		self.lines1 = []

	def __call__(self, event):
		if event.inaxes == self.ax1:
			self.process1(event)
		else:
			pass
		
		return True
		
	def process1(self, event):
		x = event.xdata
		if self.nClick1 == 0:
			self.x11 = x

			for line in self.lines1:
				try:
					self.ax1.lines.remove(line)
				except:
					pass
			self.lines1 = []
			
		else:
			self.x12 = x
			
		line, = self.ax1.plot([x,x], [-10, 10], linestyle='--', color='black')
		self.lines1.append(line)
		self.ax1.set_ylim(-10, 10)
		plt.draw()

		self.nClick1 += 1
		if self.nClick1 == 2:
			if self.x12 < self.x11:
				temp = self.x11
				self.x11 = self.x12
				self.x12 = temp
			
			toUse = numpy.where( (self.freq >= self.x11) & (self.freq <= self.x12) )[0]
			partF = self.freq[toUse]
			if self.extra is None:
				partP = self.phase[toUse]
			else:
				if self.extra.corr is None:
					partP = self.phase[toUse]
				else:
					partP = self.extra.corr[toUse]
			
			#coeff = numpy.polyfit(partF, partP, 1)
			coeff = robust.polyfit(partF, partP, 1)
			print coeff[0] / (2*numpy.pi) * 1e9
			self.delay = coeff[0] / (2*numpy.pi)
			
			line, = self.ax1.plot(freq, numpy.polyval(coeff, freq), linestyle=':', color='red')
			self.ax1.set_ylim(-10, 10)
			self.lines1.append(line)
			
			if self.extra is None:
				try:
					del self.ax2.lines[0]
					del self.ax2.lines[0]
					del self.ax2.lines[0]
					del self.ax2.lines[0]
				except:
					pass
				corr =  numpy.angle( numpy.exp(1j*self.phase) / numpy.exp(1j*numpy.polyval(coeff, freq)) )
			else:
				try:
					del self.ax2.lines[3]
				except:
					pass
				corr = numpy.angle( numpy.exp(1j*self.extra.corr) / numpy.exp(1j*numpy.polyval(coeff, freq)) )
			self.corr = corr
			self.ax2.plot(freq, self.corr)
			self.ax2.set_ylim(-10, 10)
			
			plt.draw()
			
			self.nClick1 = 0


antennas = lwa1.getAntennas()

dataDict = numpy.load('prepared-dat.npz')
refAnt = dataDict['refAnt'].item()
refX   = dataDict['refX'].item()
refY   = dataDict['refY'].item()

freq = dataDict['freq']
time = dataDict['time']
data = dataDict['data']
simPhase = dataDict['simPhase']

#
# Downselect frequency coverage
#
good = numpy.where( freq <= 80e6 )[0]
freq = freq[good]
time = time[good]
data = data[good,:]
simPhase = simPhase[good,:]

print freq.shape

#
# Apply model correction
#
data /= simPhase

#
# Select antenna to work on
#
try:
	toWork = int(sys.argv[1])
except:
	toWork = 0


#
# Plot
#
fig = plt.figure()
ax1 = fig.add_subplot(1, 2, 1)
ax2 = fig.add_subplot(1, 2, 2)
ax1.plot(freq, numpy.angle( data[:,toWork] ))
ax1.set_ylim(-10, 10)
ax1.set_title('Antenna #%i, Status %i' % (antennas[toWork].id, antennas[toWork].getStatus()))
ax2.set_title('Stand #%i, Pol. %s' % (antennas[toWork].stand.id, 'X' if antennas[toWork].pol == 0 else 'Y'))
ax2.grid(True)
fr1 = FitRegion(ax1, ax2, freq, numpy.angle( data[:,toWork] ))
fr2 = FitRegion(ax2, ax2, freq, numpy.angle( data[:,toWork] ), extra=fr1)
plt.show()

#
# Combine
#
centralFreq = numpy.median(freq)
delay = fr1.delay + fr2.delay
print antennas[refX].cable.delay(centralFreq) - delay, antennas[refX].cable.delay(centralFreq) - delay - antennas[toWork].cable.delay(centralFreq)

#
# Save
#
fh = open('add-delay.txt', 'a')
fh.write('%3i  %.6g\n' % (toWork, antennas[refX].cable.delay(centralFreq) - delay - antennas[toWork].cable.delay(centralFreq)))
fh.close()

