#!/usr/bin/env python

"""
Simulate fringes for a dipole-dipole data set using the lsl.sim.vis.build_sim_data()
function and the bright sources listed in lsl.sim.vis.srcs.
"""

# Python3 compatiability
from __future__ import print_function, division
import sys
if sys.version_info > (3,):
    xrange = range
    
import os
import sys
import numpy
import argparse
from calendar import timegm
from datetime import datetime

import aipy

from lsl.reader import drx
from lsl.common.dp import fS
from lsl.common import stations
from lsl.common.paths import DATA as dataPath
from lsl.astro import unix_to_utcjd
from lsl.sim import vis as simVis

from matplotlib import pyplot as plt

def main(args):
    filenames = args.filenames
    filenames.sort()

    times = []
    for filename in filenames:
        dataDict = numpy.load(filename)

        tStart = datetime.utcfromtimestamp(dataDict['tStart'])
        tInt = dataDict['tInt']
        try:
            srate = dataDict['srate']
        except KeyError:
            srate = 19.6e6
        
        freq1 = dataDict['freq1']
        freq2 = dataDict['freq2']
        stand1, stand2 = dataDict['stands']

        times.append( tStart)

    print("Got %i files from %s to %s (%s)" % (len(filenames), times[0].strftime("%Y/%m/%d %H:%M:%S"), times[-1].strftime("%Y/%m/%d %H:%M:%S"), (times[-1]-times[0])))

    iTimes = []
    for i in xrange(1, len(times)):
        dt = times[i] - times[i-1]
        iTimes.append(dt.days*24*3600 + dt.seconds + dt.microseconds/1e6)
    iTimes = numpy.array(iTimes)
    print(" -> Interval: %.3f +/- %.3f seconds (%.3f to %.3f seconds)" % (iTimes.mean(), iTimes.std(), iTimes.min(), iTimes.max()))
    
    print("Number of frequency channels: %i (~%.1f Hz/channel)" % (len(freq1)+1, freq1[1]-freq1[0]))

    # Build up the station
    if args.lwasv:
        site = stations.lwasv
    else:
        site = stations.lwa1
    
    rawAntennas = site.antennas
    
    antennas = []
    for ant in rawAntennas:
        if ant.stand.id == stand1 and ant.pol == 0:
            antennas.append(ant)
    for ant in rawAntennas:
        if ant.stand.id == stand2 and ant.pol == 0:
            antennas.append(ant)
    if len(antennas) != 2:
        raise RuntimeError("Can only find stand %i, %i and %i found in the NPZ files" % (antennas[0].stand.id, stand1, stand2))

    # Create the simulated array
    refJD = unix_to_utcjd(timegm(times[0].timetuple()))
    aa1 = simVis.build_sim_array(site, antennas, freq1/1e9, jd=refJD)
    aa2 = simVis.build_sim_array(site, antennas, freq2/1e9, jd=refJD)

    # Build the model times and range.
    jdList = []
    dTimes = []
    for i in xrange(len(times)):
        tNow = timegm(times[i].timetuple())
        jdNow = unix_to_utcjd(tNow)

        jdList.append(jdNow)
        dTimes.append( (times[i]-times[0]).seconds )
        
    # Actually run the simulations
    simDict1 = simVis.build_sim_data(aa1, simVis.SOURCES, jd=jdList, pols=['xx',], verbose=False)
    simDict2 = simVis.build_sim_data(aa2, simVis.SOURCES, jd=jdList, pols=['xx',], verbose=False)

    # Plot
    fig = plt.figure()
    ax1 = fig.add_subplot(2, 1, 1)
    ax2 = fig.add_subplot(2, 1, 2)

    vis1 = []
    for vis in simDict1['vis']['xx']:
        vis1.append( vis )
    vis2 = []
    for vis in simDict2['vis']['xx']:
        vis2.append( vis )
    
    vis1 = numpy.array(vis1)
    vis1 = numpy.ma.array(vis1, mask=~numpy.isfinite(vis1))
    vis2 = numpy.array(vis2)
    vis2 = numpy.ma.array(vis2, mask=~numpy.isfinite(vis2))

    data = numpy.abs(vis1)
    data = data.ravel()
    data.sort()
    vmin1 = data[int(round(0.15*len(data)))]
    vmax1 = data[int(round(0.85*len(data)))]
    print('Plot range for tuning 1:', vmin1, vmax1)
    
    data = numpy.abs(vis2)
    data = data.ravel()
    data.sort()
    vmin2 = data[int(round(0.15*len(data)))]
    vmax2 = data[int(round(0.85*len(data)))]
    print('Plot range for tuning 2:', vmin2, vmax2)

    ax1.imshow(numpy.abs(vis1), extent=(freq1[0], freq1[-1], dTimes[0], dTimes[-1]), origin='lower', 
            vmin=vmin1, vmax=vmax1)
    ax2.imshow(numpy.abs(vis2), extent=(freq1[0], freq1[-1], dTimes[0], dTimes[-1]), origin='lower', 
            vmin=vmin2, vmax=vmax2)

    ax1.axis('auto')
    ax2.axis('auto')

    fig.suptitle("%s to %s UTC" % (times[0].strftime("%Y/%m/%d %H:%M"), times[-1].strftime("%Y/%m/%d %H:%M")))
    ax1.set_xlabel('Frequency [MHz]')
    ax2.set_xlabel('Frequency [MHz]')
    ax1.set_ylabel('Elapsed Time [s]')
    ax2.set_ylabel('Elapsed Time [s]')

    plt.show()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
    description='Simulate fringes for a dipole-dipole data set using the lsl.sim.vis.build_sim_data() function and the bright sources listed in lsl.sim.vis.srcs.',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('filenames', type=str, nargs='+',
            help='Data sets (.npz)')
    parser.add_argument('-v','--lwasv', action='store_true',
            help='Station is LWA-SV')

    args = parser.parse_args()
    main(args)
