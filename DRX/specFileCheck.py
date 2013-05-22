#!/usr/bin/env python

# Description
# This script takes a DROSv2  XX/YY data and checks the power and clip levels

import os
import sys
import h5py
import numpy
import getopt
from datetime import datetime

from lsl.reader import drspec
from lsl.common import progress
from optparse import OptionParser

if __name__ == '__main__':
    parser = OptionParser(usage="inspect_data.py [file]", description="Inspect data in file")
    (options, args) = parser.parse_args()
    
    fh = open(args[0], 'rb')
    # Interogate the file to figure out what frames sizes to expect, now many 
    # frames there are, and what the transform length is
    nFrames = os.path.getsize(args[0]) / drspec.getFrameSize(fh)
    nChunks = nFrames
    LFFT = drspec.getTransformSize(fh)
    
    # Read in the first frame to figure out the DP information
    cPos = fh.tell()
    junkFrame = drspec.readFrame(fh)
    fh.seek(cPos)
    
    beam = junkFrame.parseID()
    centralFreq1 = junkFrame.getCentralFreq(1)
    centralFreq2 = junkFrame.getCentralFreq(2)
    srate = junkFrame.getSampleRate()
    tInt = junkFrame.header.nInts*LFFT/srate
    decFactor = junkFrame.header.decimation
    beginDate = datetime.utcfromtimestamp(junkFrame.getTime())
    
    products = junkFrame.getDataProducts()
    nProducts = len(products)
    
    # Report
    print "Filename: %s" % args[0]
    print "Date of First Frame: %s" % beginDate
    print "Beam: %i" % beam
    print "Sample Rate: %i Hz" % srate
    print "Tuning Frequency: %.3f Hz (1); %.3f Hz (2)" % (centralFreq1, centralFreq2)
    print "Data Products: %s" % (', '.join(products),)
    print "Frames: %i (%.3f s)" % (nFrames, nFrames*tInt)
    print "---"
    print "Transform Length: %i" % LFFT
    print "Integration: %.3f s" % tInt
    print "Decimation Factor: %.3f" % decFactor

    percentages = numpy.zeros((4,nFrames),numpy.dtype('float'))
    power = numpy.zeros((2*nProducts,nFrames),numpy.dtype('float'))
 
    if products == ['XX', 'YY'] or products == ['XX','YY','XY','YX']: 
    	for i in xrange(nFrames):
        	frame = drspec.readFrame(fh)
		if nProducts == 4:
			data = [ frame.data.XX0.mean(),  frame.data.YY0.mean(),  frame.data.XY0.mean(),  frame.data.YX0.mean(), frame.data.XX1.mean(), frame.data.YY1.mean(), frame.data.XY1.mean(), frame.data.YX1.mean() ]
		else:
        		data = [ frame.data.XX0.mean(),  frame.data.YY0.mean(),  frame.data.XX1.mean(), frame.data.YY1.mean() ]
        	for j in range(nProducts*2):
            		power[j][i] = data[j]
		for j in range(4):
            		percentages[j][i] = int(frame.data.saturations[j])/(srate*tInt)
        
        
    	print "Average Saturation Count in %:"
	if nProducts == 2:
    		print "Tuning 1 XX: %f (%f)" % (percentages[0].mean()*100,  percentages[0].std()*100)
    		print "Tuning 1 YY: %f (%f)" % (percentages[1].mean()*100,  percentages[1].std()*100)
    		print "Tuning 2 XX: %f (%f)" % (percentages[2].mean()*100,  percentages[2].std()*100)
    		print "Tuning 2 YY: %f (%f)" % (percentages[3].mean()*100,  percentages[3].std()*100)
    		print "---"
    		print "Average Power Levels:"
	if nProducts == 2:
    		print "Tuning 1 XX: %f" % (power[0].mean())
    		print "Tuning 1 YY: %f" % (power[1].mean())
    		print "Tuning 2 XX: %f" % (power[2].mean())
    		print "Tuning 2 YY: %f" % (power[3].mean())
	else:
                print "Tuning 1 XX: %f" % (power[0].mean())
                print "Tuning 1 YY: %f" % (power[1].mean())
		print "Tuning 1 XY: %f" % (power[2].mean())
                print "Tuning 1 YX: %f" % (power[3].mean())
		print "Tuning 2 XX: %f" % (power[4].mean())
                print "Tuning 2 YY: %f" % (power[5].mean())
		print "Tuning 2 XY: %f" % (power[6].mean())
                print "Tuning 2 YX: %f" % (power[7].mean())

    elif products == ['I','Q','U','V'] or products== ['I','V']:
	for i in xrange(nFrames):
		frame = drspec.readFrame(fh)
		if nProducts == 4:
			data = [ frame.data.I0.mean(), frame.data.Q0.mean(), frame.data.U0.mean(), frame.data.V0.mean(), frame.data.I1.mean(), frame.data.Q1.mean(), frame.data.U1.mean(), frame.data.V1.mean() ]
		else:
			data = [ frame.data.I0.mean(), frame.data.V0.mean(), frame.data.I1.mean(), frame.data.V1.mean() ]
		for j in range(nProducts*2):
			power[j][i] = data[j]
		for j in range(4):
			percentages[j][i] = int(frame.data.saturations[j])/(srate*tInt)

	print "Average Saturation Count in %:"
	print "Tuning 1 XX: %f (%f)" % (percentages[0].mean()*100, percentages[0].std()*100)
	print "Tuning 1 YY: %f (%f)" % (percentages[1].mean()*100, percentages[1].std()*100)
	print "Tuning 2 XX: %f (%f)" % (percentages[2].mean()*100, percentages[2].std()*100)
	print "Tuning 2 YY: %f (%f)" % (percentages[3].mean()*100, percentages[3].std()*100)
 	print "---"
	print "Average Power Levels"
	if nProducts==2:
		print "Tuning 1 I: %f" % (power[0].mean())
                print "Tuning 1 V: %f" % (power[1].mean())
                print "Tuning 2 I: %f" % (power[2].mean())
                print "Tuning 2 V: %f" % (power[3].mean())

	else:
        	print "Tuning 1 I: %f" % (power[0].mean())
        	print "Tuning 1 Q: %f" % (power[1].mean())
        	print "Tuning 1 U: %f" % (power[2].mean())
        	print "Tuning 1 V: %f" % (power[3].mean())
        	print "Tuning 2 I: %f" % (power[4].mean())
        	print "Tuning 2 Q: %f" % (power[5].mean())
        	print "Tuning 2 U: %f" % (power[6].mean())
        	print "Tuning 2 V: %f" % (power[7].mean())

    else:
        print "For this spectrometer mode not implemented"
  
    fh.close()
    
    
