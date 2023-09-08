#!/usr/bin/env python3

import os
import sys
import aipy
import ephem
import numpy
import argparse

from datetime import datetime
from scipy.optimize import leastsq

from lsl.common.stations import lwa1
from lsl.sim.vis import SOURCES as simSrcs

from matplotlib import pyplot as plt


# List of bright radio sources and pulsars in PyEphem format
_srcs = ["TauA,f|J,05:34:32.00,+22:00:52.0,1", 
         "VirA,f|J,12:30:49.40,+12:23:28.0,1",
         "CygA,f|J,19:59:28.30,+40:44:02.0,1", 
         "CasA,f|J,23:23:27.94,+58:48:42.4,1",
         "3C123,f|J,04:37:04.38,+29:40:13.8,1",
         "3C295,f|J,14:11:20.47,+52:12:09.5,1",
         "HerA,f|J,16:51:08.15,+04:59:33.3,1",
         "SgrA,f|J,17:45:40.00,-29:00:28.0,1"]


class RadioBodyBaars:
    """
    Class defining flux of a celestial source using the Baars a, b, and c parameters.

    Based on the aipy.amp.RadioBody class
    """
    def __init__(self, a, b, c, secularChange=0.0, secularEpoch=2013.0):
        """
        Flux parameters where:
        log S = a + b*log(nu/1MHz) + c*log(nu/1MHz)**2
        with the posibility of a secular evolution since a certain
        epoch (as a year).
        """
        
        self.a = a
        self.b = b
        self.c = c
        self.secularChange = secularChange
        self.secularEpoch = secularEpoch
        self.mfreq = 1e6 / 1e9	# 1 MHz
        
    def update_jys(self, afreqs, epoch=2013.0):
        """
        Update fluxes relative to the provided observer.  Must be called at 
        each time step before accessing information.
        """

        flux = 10**(self.a + self.b*numpy.log10(afreqs/self.mfreq) + self.c*numpy.log10(afreqs/self.mfreq)**2)
        flux *= (1 + self.secularChange)**(epoch-self.secularEpoch)
        
        self.jys = flux
        
    def get_jys(self):
        """
        Return the fluxes vs. freq that should be used for simulation.
        """
        
        return self.jys


class RadioFixedBodyBaars(aipy.phs.RadioFixedBody, RadioBodyBaars):
    """
    Class representing a source at fixed RA,DEC.  Adds Baars-style flux 
    information to aipy.phs.RadioFixedBody.

    Based on the aipy.amp.RadioFixedBody class
    """
    
    def __init__(self, ra, dec, name='', epoch=ephem.J2000, a=1.0, b=0.0, c=0.0, secularChange=0.0, secularEpoch=2013.0, mfreq=0.001, ionref=(0.,0.), srcshape=(0.,0.,0.), **kwargs):
        """
        ra = source's right ascension (epoch=J2000)
        dec = source's declination (epoch=J2000)
        jys = source strength in Janskies at mfreq)
        mfreq = frequency (in GHz) where source strength was measured
        index = power-law index of source emission vs. freq.
        """
        
        aipy.phs.RadioFixedBody.__init__(self, ra, dec, mfreq=mfreq, name=name, epoch=epoch, ionref=ionref, srcshape=srcshape)
        RadioBodyBaars.__init__(self, a, b, c, secularChange=secularChange, secularEpoch=secularEpoch)
        
    def compute(self, observer, afreqs=74e-3):
        epoch = datetime.strptime(str(observer.date), "%Y/%m/%d %H:%M:%S")
        epoch = epoch.year + float(epoch.strftime("%j")) / 365

        aipy.phs.RadioFixedBody.compute(self, observer)
        try:
            self.update_jys(observer.get_afreqs(), epoch=epoch)
        except AttributeError:
            self.update_jys(afreqs, epoch=epoch)


def getSources():
    """
    Return a dictionary of PyEphem sources.
    """
    
    srcs = {}
    for line in _srcs:
        src = ephem.readdb(line)
        srcs[src.name] = src
        
    return srcs


def getAIPYSources():
    """
    Return a dictionary of AIPY sources.
    
    .. note::
        This function returns a slightly different dictionary that than 
        contained in lsl.sim.vis.srcs.  First, this dictionary has source 
        names that are consistent with those returned by the getSources()
        function.  Second, this list contains 3C123.
    """
    
    newSrcs = {}
    for name,src in simSrcs.items():
        if name == 'Sun':
            newSrcs[name] = src
        elif name == 'Jupiter':
            newSrcs[name] = src
        elif name == 'crab':
            newSrcs['TauA'] = src
        else:
            newSrcs['%sA' % name.capitalize()] = src
    newSrcs['3C123'] = aipy.amp.RadioFixedBody('4:37:04.38', '+29:40:13.8', jys=206.0, index=-0.70, mfreq=0.178)
    
    # Modify CygA, TauA, VirA, and CasA for the Baars et al. (1977) fluxes.  
    # For CasA, include the Helmboldt & Kassim secular decrease of 0.84%/yr
        # -> (Helmboldt & Kassim 2009, AJ, 138, 838)
    newSrcs['CygA'] = RadioFixedBodyBaars('19:59:28.30', '+40:44:02.0', a=4.695, b=0.085, c=-0.178)
    newSrcs['TauA'] = RadioFixedBodyBaars('05:34:32.00', '+22:00:52.0', a=3.915, b=-0.299)
    newSrcs['VirA'] = RadioFixedBodyBaars('12:30:49.40', '+12:23:28.0', a=5.023, b=-0.856)
    newSrcs['CasA'] = RadioFixedBodyBaars('23:23:27.94', '+58:48:42.4', a=5.625, b=-0.634, c=-0.023,
                        secularChange=-0.0084, secularEpoch=1965.0)
    newSrcs['3C295'] = RadioFixedBodyBaars('14:11:20.47', '+52:12:09.5', a=1.485, b=0.759, c=-0.255)
    
    return newSrcs


def _driftscanFunction(p, x):
    height = p[0]
    center = p[1]
    width  = p[2]
    offset = p[3]
    try:
        slope  = p[4]
    except IndexError:
        slope = 0.0
    
    y = height*numpy.exp(-4*numpy.log(2)*(x - center)**2/width**2 ) + slope*x + offset
    return y


def _driftscanErrorFunction(p, x, y):
    yFit = _driftscanFunction(p, x)
    return y - yFit 


def fitDriftscan(t, power, includeLinear=False):
    """
    Given an array of times and and array of total power from a drift scan, 
    fit the drift scan with a Gaussian to estimate the RA error, SEFD, and
    FWHM.  Return the results as a four-element tuple of:
    1) time of peak, 
    2) SEFD metric (1/(P1/P0 - 1), 
    3) FWHM in seconds
    4) best-fit values
    """
    
    gp = [power.max()-power.min(), t.mean(), 1000, power.min()]
    if includeLinear:
        gp.append(0.0)
    gp, status = leastsq(_driftscanErrorFunction, gp, (t, power))
    
    tPeak = gp[1]
    sefdMetric = gp[3]/gp[0]
    fwhm = abs(gp[2])
    
    fit = _driftscanFunction(gp, t)
    
    if includeLinear:
        slope = gp[4]
        return tPeak, sefdMetric, fwhm, slope, fit
    else:
        return tPeak, sefdMetric, fwhm, fit


def _decFunction(p, x, fwhm):
    height = p[0]
    center = p[1]
    offset = p[2]
    width = fwhm
    
    y = height*numpy.exp(-4*numpy.log(2)*(x - center)**2/width**2 ) + offset
    return y


def _decErrorFunction(p, x, y, fwhm):
    yFit = _decFunction(p, x, fwhm)
    return y - yFit 


def fitDecOffset(decs, powers, fwhm=2.0):
    """
    Given an array of declination offsets from the source and the peak 
    power seen from each driftscan, estimate the declination pointing 
    error.
    
    .. note::
        To reduce the number of samples needed for the fit, the FWHM is
        specified during the fitting.  The default value used is two 
        degrees.
    """
    
    gp = [powers.max()-powers.min(), 0.0, powers.min()]
    gp, status = leastsq(_decErrorFunction, gp, (decs, powers, fwhm))
    
    return gp[1]


def main(args):
    # Parse the command line
    filename = args.filename
    
    # Load in the data to see what we should do
    dataDict = numpy.load(filename)
    srcName = dataDict['source'].item()
    freq = dataDict['freq']
    unx = dataDict['unx']
    lst = dataDict['lst']
    pwrX = dataDict['pwrX']
    pwrY = dataDict['pwrY']
    dataDict.close()
    
    # Form Stokes I out of X and Y
    pwrI = pwrX + pwrY
    
    # Read in each of the data sets
    data = {}
    pointing = {}
    tStartPlot = 1e12
    for i,name in enumerate(('south', srcName, 'north')):
        pointing[name] = name
        
        if unx[0] < tStartPlot:
            tStartPlot = unx[0]
            
        data[name] = {'t':unx, 'f1':freq, 'I1': pwrI[:,i,:]}
        
    # Get LWA-1
    observer = lwa1.get_observer()
    
    # Load in the sources and find the right one
    srcs = getSources()
    simSrcs = getAIPYSources()
    
    toUse = None
    for srcName in data.keys():
        for src in srcs.keys():
            if src.lower() == srcName.lower():
                toUse = src
                observer.date = datetime.utcfromtimestamp( data[srcName]['t'][0] )
                break;
    if toUse is None:
        raise RuntimeError("Unknown source in input files")
        
    toUseAIPY = srcs[toUse].name
    try:
        simSrcs[toUseAIPY]
    except KeyError:
        toUseAIPY = None
        print("Warning: Cannot find flux for this target")
        
    # Find out when the source should have transitted the beam
    tTransit = 0.0
    zenithAngle = ephem.degrees('180:00:00')
    for name in pointing.keys():
        if name.lower() != srcs[toUse].name.lower():
            continue
            
        observer.next_transit(srcs[toUse])
        az = srcs[toUse].az
        el = srcs[toUse].alt
        
        bestT = 0.0
        bestV = 1e6
        for t in data[name]['t']:
            observer.date = datetime.utcfromtimestamp(t).strftime("%Y/%m/%d %H:%M:%S")
            srcs[toUse].compute(observer)
            
            sep = ephem.separation((srcs[toUse].az, srcs[toUse].alt), (az,el))
            if sep < bestV:
                bestT = t
                bestV = sep
    tTransit = bestT
    observer.date = datetime.utcfromtimestamp(tTransit).strftime("%Y/%m/%d %H:%M:%S")
    zenithAngle = ephem.degrees(ephem.degrees('90:00:00') - el)
    
    # Plot
    fig = plt.figure()
    ax1 = fig.gca()
    
    raOffsets1 = {}
    decPowers1 = {}
    fwhmEstimates1 = {}
    sefdEstimate1 = None
    
    raOffsets2 = {}
    decPowers2 = {}
    fwhmEstimates2 = {}
    sefdEstimate2 = None
    for name in data.keys():
        t = data[name]['t']
        f1 = data[name]['f1']
        I1 = data[name]['I1']
        
        # Select data that was actually recorded
        good = numpy.where( (t > 0) & (I1[:,10] > 0) )[0][:-1]
        t = t[good]
        I1 = I1[good]
        
        # Sum over the interesting part of the band
        toUseSpec = numpy.where( (f1 > 60e6) & (f1 < 80e6) )[0]
        I1 = I1[:,toUseSpec].sum(axis=1)
        
        # Convert the scales to unit flux
        I1 /= I1.max()
        
        # Fit a Gaussian to the power to find the transit time and beam width
        includeLinear = False
        if includeLinear:
            obsTransit1, sefdMetric1, obsFWHM1, obsSlope1, obsFit1 = fitDriftscan(t, I1, includeLinear=True)
            
            linear1 = obsSlope1 * t
            linear1 -= linear1[:10].mean()
            
            I1 -= linear1
            obsFit1 -= linear1
            
        obsTransit1, sefdMetric1, obsFWHM1, obsFit1 = fitDriftscan(t, I1, includeLinear=False)
        
        # Save the results
        diff1 = obsTransit1 - tTransit
        raOffsets1[name] = diff1
        decPowers1[name] = obsFit1.max() - obsFit1.min()
        fwhmEstimates1[name] = obsFWHM1/3600.*15.*numpy.cos(srcs[toUse]._dec)
        
        # Report
        print('Target: %s' % name)
        print('  Tuning 1 @ %.2f MHz' % (f1.mean()/1e6,))
        print('    FWHM: %.2f s (%.2f deg)' % (obsFWHM1, obsFWHM1/3600.*15.*numpy.cos(srcs[toUse]._dec)))
        print('    Observed Transit: %s' % datetime.utcfromtimestamp(obsTransit1))
        print('    Expected Transit: %s' % datetime.utcfromtimestamp(tTransit))
        print('    -> Difference: %.2f s' % diff1)
        if toUseAIPY is None:
            print('    1/(P1/P0 - 1): %.3f' % sefdMetric1)
        else:
            simSrcs[toUseAIPY].compute(observer, afreqs=f1.mean()/1e9)
            srcFlux = simSrcs[toUseAIPY].jys
            sefd = srcFlux*sefdMetric1 / 1e3
            print('    S / (P1/P0 - 1): %.3f kJy' % sefd)
            if name == srcs[toUse].name:
                sefdEstimate1 = sefd*1e3
                
        # Plot
        ax1.plot(t-tStartPlot, I1, label="%s" % name)
        ax1.plot(t-tStartPlot, obsFit1, linestyle=':')
        
    ylim1 = ax1.get_ylim()
    ax1.vlines(tTransit-tStartPlot, *ylim1, linestyle='--', label='Expected Transit')
    ax1.set_ylim(ylim1)
    ax1.legend(loc=0)
    ax1.set_title('%.2f MHz' % (f1.mean()/1e6,))
    ax1.set_xlabel('Elapsed Time [s]')
    ax1.set_ylabel('Power [arb.]')
    
    # Compute the dec. offset
    dataSet1 = (f1.mean(), raOffsets1, decPowers1, fwhmEstimates1, sefdEstimate1)
    
    sys.stderr.write("Source YYYY/MM/DD HH:MM:SS MHz    Z          errRA      errDec      SEFD      FWHM\n")
    for f,raOffsets,decPowers,fwhmEstimates,sefdEstimate in (dataSet1,):
        do = []
        dp = []
        bestOffset = None
        bestPower = -1e6
        bestFWHM = None
        for name in decPowers.keys():
            offset = raOffsets[name]
            power = decPowers[name]
            fwhm = fwhmEstimates[name]
        
            if power > bestPower:
                bestOffset = offset
                bestPower = power
                bestFWHM = fwhm
            
            if name.find('north') != -1:
                do.append(1.0)
            elif name.find('south') != -1:
                do.append(-1.0)
            else:
                do.append(0.0)
            dp.append(power)
        
        do = numpy.array(do)
        dp = numpy.array(dp)
        order = numpy.argsort(do)
        do = do[order]
        dp = dp[order]
    
        try:
            decOffset = fitDecOffset(do, dp, fwhm=bestFWHM)
        except TypeError:
            decOffset = -99
        
        raOffset = ephem.hours('00:00:%f' % bestOffset)
        decOffset = ephem.degrees('%f' % decOffset)
        fwhmEstimate = ephem.degrees('%f' % bestFWHM)
        print("%-6s %-19s %6.3f %-10s %-10s %-10s %10.3f %-10s" % (srcs[toUse].name, datetime.utcfromtimestamp(tTransit).strftime("%Y/%m/%d %H:%M:%S"), f/1e6, zenithAngle, raOffset, decOffset, sefdEstimate, fwhmEstimate))
        
    plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="fit a Gaussian to the results of estimateSEFD.py to determine the SEFD and pointing error",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
    parser.add_argument('filename', type=str,
                        help='filename to process')
    args = parser.parse_args()
    main(args)
    
