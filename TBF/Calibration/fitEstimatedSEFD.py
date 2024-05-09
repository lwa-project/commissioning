#!/usr/bin/env python3

import os
import sys
import aipy
import ephem
import numpy as np
import argparse

from datetime import datetime
from scipy.optimize import leastsq

from lsl.common.stations import parse_ssmif

from analysis import getSources, getAIPYSources, fitDriftscan, fitDecOffset

from matplotlib import pyplot as plt


def func(p, x):
    if len(p) == 4:
        height = p[0]
        center = p[1]
        width  = p[2]
        offset = p[3]
    else:
        height = p[0]
        center = p[1]
        width  = 2.0
        offset = p[2]
    
    y = height*np.exp(-4*np.log(2)*(x - center)**2/width**2 ) + offset
    return y


def err(p, x, y):
    yFit = func(p, x)
    return y - yFit 


def main(args):
    # Parse the command line
    filename = args.filename
    finalResults = []
    
    # Load in the data to see what we should do
    dataDict = np.load(filename)
    srcName = dataDict['source'].item()
    freq = dataDict['freq'].item()
    unx = dataDict['unx']
    lst = dataDict['lst']
    pnts = dataDict['pnts']
    pwrX = dataDict['pwrX']
    pwrY = dataDict['pwrY']
    dataDict.close()
    
    # Form Stokes I out of X and Y
    pwrI = pwrY# + pwrY
    
    # Get an observer
    _, ssmif, _ = os.path.basename(filename).split('-', 2)
    sta = parse_ssmif(ssmif+'.txt')
    sta_name = sta.name
    observer = sta.get_observer()
    
    # Load in the sources and find the right one
    srcs = getSources()
    simSrcs = getAIPYSources()
    
    toUse = None
    print(srcName)
    for src in srcs.keys():
        if src.lower() == srcName.lower():
            toUse = src
            observer.date = datetime.utcfromtimestamp( unx[unx.size//2] )
            break
    if toUse is None:
        raise RuntimeError("Unknown source in input files")
    raCtr = srcs[toUse]._ra*12/np.pi
    decCtr = srcs[toUse]._dec*180/np.pi
    
    toUseAIPY = srcs[toUse].name
    try:
        simSrcs[toUseAIPY]
    except KeyError:
        toUseAIPY = None
        print("Warning: Cannot find flux for this target")
        
    # Find out when the source should have transitted the beam
    tTransit = 0.0
    zenithAngle = ephem.degrees('180:00:00')
    observer.date = datetime.utcfromtimestamp(unx[0]).strftime("%Y/%m/%d %H:%M:%S")
    srcs[toUse].compute(observer)
    az = srcs[toUse].az
    el = srcs[toUse].alt
    rpos = f"azimuth: {az*180/np.pi:.1f}, elevation: {el*180/np.pi:1f}"
    
    bestT = 0.0
    bestV = 1e6
    for v in unx:
        observer.date = datetime.utcfromtimestamp(v).strftime("%Y/%m/%d %H:%M:%S")
        srcs[toUse].compute(observer)
        
        sep = ephem.separation((srcs[toUse].az, srcs[toUse].alt), (az,el))
        if sep < bestV:
            bestT = v
            bestV = sep
    tTransit = bestT
    observer.date = datetime.utcfromtimestamp(tTransit).strftime("%Y/%m/%d %H:%M:%S")
    zenithAngle = ephem.degrees(ephem.degrees('90:00:00') - el)
    
    # Convert the scales to unit flux
    pwrI /= pwrI.max()
    
    # Break the data into two pieces:
    #  raCut - Points that are part of the RA cut (second half)
    #  decCut Points that are part of the Dec cut (first half)
    raCut  = [i for i,s in enumerate(pnts) if (i > len(pnts)//2)]
    decCut = [i for i,s in enumerate(pnts) if (i < len(pnts)//2)] 
    
    # Pull out the mean data value for each step that has been observed
    m, ra, dec, pwr1 = [], [], [], []
    for i in range(len(pnts)):
        m.append( unx[0] )
        ra.append( pnts[i][0]*12/np.pi )
        dec.append( pnts[i][1]*180/np.pi )
        pwr1.append( pwrI[0,i] )
    m, ra, dec, pwr1 = np.array(m), np.array(ra), np.array(dec), np.array(pwr1)
    
    # Weed out any bad points in the RA and dec cuts
    ## RA
    while True:
        for i,p in enumerate(raCut):
            if not np.isfinite(pwr1[p]):
                del raCut[i]
                break
        break
    ## Dec
    while True:
        for i,p in enumerate(decCut):
            if not np.isfinite(pwr1[p]):
                del decCut[i]
                break
        break
        
    # Plots and analysis
    fig = plt.figure()
    fig.suptitle('Source: %s @ %s\n%s' % (srcName, sta_name, rpos))
    ax11 = fig.add_subplot(2, 1, 1)
    ax12 = fig.add_subplot(2, 1, 2)
    for i,(ax1,ax2),f,pwr in zip((1,), ((ax11,ax12),), (freq,), (pwr1,)):
        if i == 2 and tuning1 is tuning2:
            continue
        print("Tuning %i @ %.3f MHz" % (i, f/1e6))
        
        ## Dec
        x = dec[decCut]
        xPrime = np.linspace(x.min(), x.max(), 101)
        y = pwr[decCut]
        p0 = (y.max()-y.min(), x.mean(), 2.0, y.min())
        p, status = leastsq(err, p0, args=(x, y))
        decOffset = ephem.degrees(str(p[1] - decCtr))
        fwhmD = ephem.degrees(str(p[2]))
        sefdMetricD = p[3] / p[0]
        print("  Dec")
        print("    FWHM Estimate: %s" % fwhmD)
        print("    Pointing Error: %s" % decOffset)
        if toUseAIPY is None:
            print("    1/(P1/P0 - 1): %.3f" % sefdMetricD)
            if srcName == srcs[toUse].name:
                    sefdEstimateD = np.nan
        else:
            try:
                simSrcs[toUseAIPY].compute(observer, afreqs=f/1e9)
                srcFlux = simSrcs[toUseAIPY].jys
            except TypeError:
                f0, index, Flux0 = simSrcs[toUseAIPY].mfreq, simSrcs[toUseAIPY].index, simSrcs[toUseAIPY]._jys
                srcFlux = Flux0 * (f/1e9 / f0)**index
            sefd = srcFlux*sefdMetricD / 1e3
            print("    S / (P1/P0 - 1): %.3f kJy" % sefd)
            if srcName == srcs[toUse].name:
                    sefdEstimateD = sefd*1e3
                    
        ax = ax1
        ax.plot(x, y, linestyle='', marker='+', label='Data')
        ax.plot(xPrime, func(p, xPrime), linestyle='-', label='Fit')
        ax.vlines(decCtr, *ax.get_ylim(), linestyle=':')
        ax.legend(loc=0)
        ax.set_xlabel('Dec. [$^\\circ$]')
        ax.set_ylabel('Power [arb., corr.]')
        
        ## RA
        x = ra[raCut]
        xPrime = np.linspace(x.min(), x.max(), 101)
        y = pwr[raCut]
        p0 = (y.max()-y.min(), x.mean(), 2.0/15.0, y.min())
        p, status = leastsq(err, p0, args=(x, y))
        raOffset = ephem.hours(str(p[1] - raCtr))
        fwhmR = ephem.degrees(str(p[2]*15 * np.cos(decCtr*np.pi/180.0)))
        sefdMetricR = p[3] / p[0]
        print("  RA")
        print("    FWHM Estimate: %s" % fwhmR)
        print("    Pointing Error: %s" % raOffset)
        if toUseAIPY is None:
            print("    1/(P1/P0 - 1): %.3f" % sefdMetricR)
            if srcName == srcs[toUse].name:
                    sefdEstimateR = np.nan
        else:
            try:
                simSrcs[toUseAIPY].compute(observer, afreqs=f/1e9)
                srcFlux = simSrcs[toUseAIPY].jys
            except TypeError:
                f0, index, Flux0 = simSrcs[toUseAIPY].mfreq, simSrcs[toUseAIPY].index, simSrcs[toUseAIPY]._jys
                srcFlux = Flux0 * (f/1e9 / f0)**index
            sefd = srcFlux*sefdMetricR / 1e3
            print("    S / (P1/P0 - 1): %.3f kJy" % sefd)
            if srcName == srcs[toUse].name:
                    sefdEstimateR = sefd*1e3
        
        ax = ax2
        ax.plot(x, y, linestyle='', marker='+', label='Data')
        ax.plot(xPrime, func(p, xPrime), linestyle='-', label='Fit')
        ax.vlines(raCtr, *ax.get_ylim(), linestyle=':')
        ax.legend(loc=0)
        ax.set_xlabel('RA [$^h$]')
        ax.set_ylabel('Power [arb., corr.]')
        
        # Save
        fwhmEstimate = ephem.degrees((fwhmD + fwhmR) / 2.0)
        sefdEstimate = (sefdEstimateD + sefdEstimateR) / 2.0
        finalResults.append( "%-6s %-19s %6.3f %-10s %-10s %-10s %10.3f %-10s" % \
                            (srcs[toUse].name, datetime.utcfromtimestamp(tTransit).strftime("%Y/%m/%d %H:%M:%S"), f/1e6, zenithAngle, raOffset, decOffset, sefdEstimate, fwhmEstimate) )
        
    plt.show()
    
    # Final report
    sys.stderr.write("Source YYYY/MM/DD HH:MM:SS MHz    Z          errRA      errDec      SEFD      FWHM\n")
    for line in finalResults:
        print(line)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="fit a Gaussian to the results of estimateSEFD.py to determine the SEFD and pointing error",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
    parser.add_argument('filename', type=str,
                        help='filename to process')
    args = parser.parse_args()
    main(args)
    
