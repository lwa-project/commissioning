#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
plotShelterTemp.py - Script to read in the shelter.txt file and plot up the 
shelter temperature (in F) as a function of time.

$Rev$
$LastChangedBy$
$LastChangedDate$
"""

import os
import sys
import numpy
import pytz
from datetime import datetime

from matplotlib import pyplot as plt
from matplotlib import pyplot as plt
from matplotlib.dates import *
from matplotlib.ticker import *

MST7MDT = pytz.timezone('US/Mountain')


def main(args):
	if len(args) < 1:
		print 'Need a filename to plot.'
		sys.exit(1)
	filename = args[0]
	
	# Read in the data
	fh = open(filename)
	data = []
	for line in fh:
		line = line.replace('\n', '')
		fields = line.split()
		fields = [float(f) for f in fields]
		data.append( fields )

	# Split out the time and interperate it
	data = numpy.array(data)
	dates = [MST7MDT.localize(datetime.fromtimestamp(t)) for t in data[:,0]]
	print 'File spans %s to %s with %i measurements' % (dates[0], dates[-1], len(dates))
	
	# Plot
	fig = plt.figure()
	ax1 = fig.add_subplot(1, 1, 1)
	ax1.plot_date(dates, data[:,1], fmt='-', tz=MST7MDT)

	# Label and format dates
	ax1.set_title('Shelter Temperature')
	ax1.set_xlabel('Time')
	ax1.set_ylabel('Temperature [$^\circ$F]')
	fig.autofmt_xdate()

	# Show
	plt.show()


if __name__ == "__main__":
	main(sys.argv[1:])
