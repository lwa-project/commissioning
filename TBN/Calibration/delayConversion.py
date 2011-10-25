#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import numpy
from matplotlib import pyplot as plt


try:
	filename = sys.argv[1]
except:
	filename = 'add-delay.txt'

fh = open(filename)

delays = {}
for i in xrange(520):
	delays[i] = []

for line in fh:
	line = line.replace('\n', '')

	parts = line.split(None)
	ant = parts[0]
	delay = parts[1]
	ant = int(ant)
	delay = float(delay)
	
	if len(parts) == 3:
		err = float(parts[2])
		if err > 1.7e-9:
			continue
	
	delays[ant].append( delay )

fh.close()


delays2 = numpy.zeros(520)
delays3 = numpy.zeros(520) - 1000.0
for i in xrange(520):
	part = numpy.array(delays[i])
	if part.size == 0:
		print 'skipping', i
		continue
	delays2[i] = numpy.median(part)
	delays3[i] = numpy.median(part)
	print i, numpy.median(part), part.std()

suspect = numpy.where( numpy.abs(delays3) > 30e-9 )[0]
for i in suspect:
	print i, delays3[i]*1e9

fig = plt.figure()
ax = fig.gca()
good = numpy.where( delays3 != -1000 )
ax.hist(delays3[good]*1e9, bins=50)
ax.set_xlabel('Residual Delay [ns]')
ax.set_ylabel('Number of Antennas')
plt.show()


numpy.savez('add-delay.npz', delayDiffs=delays2)
