#!/usr/bin/env python3

"""
Given a DRX HDF5 waterfall file, plot it in an interactive way.
Tkinter version converted from wxPython.
"""

import os
import sys
import h5py
import math
import time
import numpy
import ephem
import argparse
import subprocess
import webbrowser
from datetime import datetime
from multiprocessing import Pool
from scipy.interpolate import interp1d
from scipy.stats import scoreatpercentile as percentile, skew, kurtosis
from scipy.signal import savgol_filter as savitzky_golay

from astropy.time import Time as AstroTime

import lsl
from lsl.common import dp
from lsl.common import stations
from lsl.reader.drx import FILTER_CODES
from lsl.misc.mathutils import to_dB, from_dB
from lsl.statistics import robust
from lsl.statistics.kurtosis import spectral_power, std as skStd
from lsl.misc import parser as aph

import tkinter as tk
from tkinter import ttk, messagebox, filedialog

import matplotlib
matplotlib.use('TkAgg')
matplotlib.interactive(True)

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.colors import Normalize
from matplotlib.collections import LineCollection
from matplotlib import cm
from matplotlib.figure import Figure

# Alias for min/max/mean/median that will work with NaNs
try:
    nan_min = numpy.nanmin
    nan_max = numpy.nanmax
except AttributeError:
    nan_min = numpy.min
    nan_max = numpy.max
try:
    nan_mean = numpy.nanmean
except AttributeError:
    nan_mean = numpy.mean
try:
    nan_median = numpy.nanmedian
except AttributeError:
    nan_median = numpy.median
try:
    nan_percentile = numpy.nanpercentile
except AttributeError:
    nan_percentile = numpy.percentile

try:
    import bottleneck
    try:
        nan_min = bottleneck.nanmin
        nan_max = bottleneck.nanmax
    except AttributeError:
        pass
    try:
        nan_mean = bottleneck.nanmean
    except AttributeError:
        pass
    try:
        nan_median = bottleneck.nanmedian
    except AttributeError:
        pass
except ImportError:
    print("WARNING: bottleneck not found, consider installing it with 'pip'")


__version__ = "0.3"
__author__ = "Jayce Dowell"


def findLimits(data, usedB=True):
    """
    Tiny function to speed up the computing of the data range for the colorbar.
    Returns a two-element list of the lowest and highest values.
    """

    dMin = nan_min(data)
    dMax = nan_max(data)

    if usedB:
        dMin = to_dB(dMin)
    if not numpy.isfinite(dMin):
        dMin = 0

    if usedB:
        dMax = to_dB(dMax)
    if not numpy.isfinite(dMax):
        dMax = dMin + 1

    return [dMin, dMax]


def bestFreqUnits(freq):
    """Given a numpy array of frequencies in Hz, return a new array with the
    frequencies in the best units possible (kHz, MHz, etc.)."""

    # Figure out how large the data are
    try:
        scale = int(math.log10(freq.max()))
    except AttributeError:
        scale = int(math.log10(freq))
    if scale >= 9:
        divis = 1e9
        units = 'GHz'
    elif scale >= 6:
        divis = 1e6
        units = 'MHz'
    elif scale >= 3:
        divis = 1e3
        units = 'kHz'
    else:
        divis = 1
        units = 'Hz'

    # Convert the frequency
    newFreq = freq / divis

    # Return units and freq
    return (newFreq, units)


class LogNorm(Normalize):
    """
    Normalize a given value to the 0-1 range on a log scale
    """

    def __call__(self, value, clip=None):
        value = numpy.ma.asarray(value)
        mask = numpy.ma.getmaskarray(value)
        value = value.filled(self.vmax+1)
        if clip:
            numpy.clip(value, self.vmin, self.vmax)

        output = (value - self.vmin) / (self.vmax - self.vmin)
        output *= 9
        output += 1
        output = numpy.log10(output)

        output = numpy.ma.array(output, mask=mask)
        if output.shape == () and not mask:
            output = int(output)  # assume python scalar
        return output


class SqrtNorm(Normalize):
    """
    Normalize a given value to the 0-1 range on a square root scale
    """

    def __call__(self, value, clip=None):
        value = numpy.ma.asarray(value)
        mask = numpy.ma.getmaskarray(value)
        value = value.filled(self.vmax+1)
        if clip:
            numpy.clip(value, self.vmin, self.vmax)

        output = (value - self.vmin) / (self.vmax - self.vmin)
        output = numpy.sqrt(output)

        output = numpy.ma.array(output, mask=mask)
        if output.shape == () and not mask:
            output = int(output)  # assume python scalar
        return output


class SqrdNorm(Normalize):
    """
    Normalize a given value to the 0-1 range on a squared scale
    """

    def __call__(self, value, clip=None):
        value = numpy.ma.asarray(value)
        mask = numpy.ma.getmaskarray(value)
        value = value.filled(self.vmax+1)
        if clip:
            numpy.clip(value, self.vmin, self.vmax)

        output = (value - self.vmin) / (self.vmax - self.vmin)
        output = output**2

        output = numpy.ma.array(output, mask=mask)
        if output.shape == () and not mask:
            output = int(output)  # assume python scalar
        return output


class AsinhNorm(Normalize):
    """
    Normalize a given value to the 0-1 range on an inverse hyperbolic sine scale
    """

    def __call__(self, value, clip=None):
        value = numpy.ma.asarray(value)
        mask = numpy.ma.getmaskarray(value)
        value = value.filled(self.vmax+1)
        if clip:
            numpy.clip(value, self.vmin, self.vmax)

        output = (value - self.vmin) / (self.vmax - self.vmin)
        output = numpy.arcsinh(output*10.0/3.0) / numpy.arcsinh(10.0/3.0)

        output = numpy.ma.array(output, mask=mask)
        if output.shape == () and not mask:
            output = int(output)  # assume python scalar
        return output


class SinhNorm(Normalize):
    """
    Normalize a given value to the 0-1 range on an hyperbolic sine scale
    """

    def __call__(self, value, clip=None):
        value = numpy.ma.asarray(value)
        mask = numpy.ma.getmaskarray(value)
        value = value.filled(self.vmax+1)
        if clip:
            numpy.clip(value, self.vmin, self.vmax)

        output = (value - self.vmin) / (self.vmax - self.vmin)
        output = numpy.sinh(output*10.0/3.0) / numpy.sinh(10.0/3.0)

        output = numpy.ma.array(output, mask=mask)
        if output.shape == () and not mask:
            output = int(output)  # assume python scalar
        return output


class HistEqNorm(Normalize):
    """
    Normalize a given value to the 0-1 range using histogram equalization
    """

    def __call__(self, value, clip=None):
        value = numpy.ma.asarray(value)
        mask = numpy.ma.getmaskarray(value)
        value = value.filled(self.vmax+1)
        if clip:
            numpy.clip(value, self.vmin, self.vmax)

        hist, bins = numpy.histogram(value, bins=256)
        hist = numpy.insert(hist, 0, 0)
        hist = hist.cumsum() / float(hist.sum())
        histeq = interp1d(bins, hist, bounds_error=False, fill_value=0.0)
        output = histeq(value)

        output = numpy.ma.array(output, mask=mask)
        if output.shape == () and not mask:
            output = int(output)  # assume python scalar
        return output


class Waterfall_GUI(object):
    def __init__(self, frame, freq=None, spec=None, tInt=None, bandpassType='data', arxFilter='split'):
        self.frame = frame
        self.press = None

        self.filename = ''
        self.index = 0
        self.filenames = None

        self.bandpass = False
        self.freq1 = freq
        self.freq2 = freq
        self.spec = spec
        self.tInt = tInt
        self.linear = False
        self.data_products = ['XX', 'YY']
        self.bandpassType = bandpassType
        self.arxFilter = arxFilter

        self.ax1a = None
        self.ax1b = None
        self.ax1c = None
        self.ax2 = None
        self.cmap = cm.get_cmap('jet')
        self.norm = Normalize

        self.spectrumClick = None
        self._keyPressCache = []

        self.bandpassCut = 0.85
        self.driftOrder  = 5
        self.driftCut    = 4
        self.kurtosisSec = 1
        self.kurtosisCut = 4

        self.frame.edited = False

    def loadData(self, filename, obsID=1):
        """
        Load in data from an NPZ file.
        """

        print("Loading file '%s'" % os.path.split(filename)[1])
        tStart = time.time()

        # Save the filename
        self.filename = filename
        self.obsID = obsID
        h = h5py.File(self.filename, 'r')
        obs = h.get('Observation%i' % obsID, None)
        if obs is None:
            raise RuntimeError("No observation #%i in file '%s'" % (obsID, os.path.basename(self.filename)))

        # Load the Data
        print(" %6.3f s - Extracting data for observation #%i" % (time.time() - tStart, obsID))
        self._bpCache = {}
        try:
            arxFilterCode = obs.attrs['ARX_Filter']

            if arxFilterCode == 0:
                self.arxFilter = 'split'
            elif arxFilterCode == 1:
                self.arxFilter = 'full'
            elif arxFilterCode == 2:
                self.arxFilter = 'reduced'

        except KeyError:
            pass
        self.beam  = obs.attrs['Beam']
        self.srate = obs.attrs['sampleRate']
        self.tInt  = obs.attrs['tInt']
        self.time  = numpy.zeros(obs['time'].shape, dtype=obs['time'].dtype)
        obs['time'].read_direct(self.time)
        try:
            fmt = obs['time'].attrs['format']
            scl = obs['time'].attrs['scale']
            try:
                fmt = fmt.decode()
                scl = scl.decode()
            except AttributeError:
                pass

            if fmt != 'unix' or scl != 'utc':
                self.time = [AstroTime(*t, format=fmt, scale=scl) for t in self.time]
                self.time = [t.utc.unix for t in self.time]
                self.time = numpy.array(self.time)

            else:
                self.time = self.time["int"] + self.time["frac"]

        except (KeyError, ValueError):
            pass

        tuning1 = obs.get('Tuning1', None)
        tuning2 = obs.get('Tuning2', None)
        if tuning2 is None:
            tuning2 = tuning1

        data_products = list(tuning1)
        mapper = {'XX': 0, 'I': 0, 'XY_real': 1, 'Q': 1, 'XY_imag': 2, 'U': 2, 'YY': 3, 'V': 3}
        for exclude in ('freq', 'Mask', 'Saturation', 'SpectralKurtosis'):
            try:
                ind = data_products.index(exclude)
                del data_products[ind]
            except ValueError:
                pass
        if data_products[0][0] in ('X', 'Y'):
            self.linear = True
            if 'XY_real' in data_products or 'XY_imag' in data_products:
                self.usedB = False
            else:
                self.usedB = True
        else:
            self.linear = False
            self.usedB = False

        # Figure out the selection process
        self.iOffset = int(round(self.frame.offset / self.tInt))
        if self.frame.duration < 0:
            self.iDuration = self.time.size - self.iOffset
        else:
            self.iDuration = int(round(self.frame.duration / self.tInt))
        ## Make sure we don't fall off the end of the file
        if self.iOffset + self.iDuration > tuning1[data_products[0]].shape[0]:
            self.iDuration = tuning1[data_products[0]].shape[0] - self.iOffset
        ## Make sure all samples have a valid time
        while self.time[self.iOffset] <= 0:
            self.iOffset += 1
        while self.time[self.iOffset+self.iDuration-1] <= 0:
            self.iDuration -= 1
        selection = numpy.s_[self.iOffset:self.iOffset+self.iDuration, :]

        if self.iOffset != 0:
            print("            -> Offsetting %i integrations (%.3f s)" % (self.iOffset, self.iOffset*self.tInt))
        print("            -> Displaying %i integrations (%.3f s)" % (self.iDuration, self.iDuration*self.tInt))

        self.time = self.time[self.iOffset:self.iOffset+self.iDuration]
        self.time -= self.time[0]

        self.freq1 = numpy.zeros(tuning1['freq'].shape, dtype=tuning1['freq'].dtype)
        tuning1['freq'].read_direct(self.freq1)
        self.freq2 = numpy.zeros(tuning2['freq'].shape, dtype=tuning2['freq'].dtype)
        tuning2['freq'].read_direct(self.freq2)

        self.spec = numpy.empty((8, self.iDuration, self.freq1.size), dtype=numpy.float32)

        for p in data_products:
            ## Tuning 1
            ind = 4*(1-1) + mapper[p]
            part = numpy.empty((self.iDuration, self.freq1.size), dtype=tuning1[p].dtype)
            tuning1[p].read_direct(part, selection)
            self.spec[ind,:,:] = part.astype(numpy.float32)

            ## Tuning 2
            ind = 4*(2-1) + mapper[p]
            part = numpy.empty((self.iDuration, self.freq2.size), dtype=tuning2[p].dtype)
            tuning2[p].read_direct(part, selection)
            self.spec[ind,:,:] = part.astype(numpy.float32)

            del part

        self.sats = numpy.empty((self.iDuration, 4), dtype=numpy.float32)
        try:
            part = numpy.empty((self.iDuration, 2), dtype=tuning1['Saturation'].dtype)
            tuning1['Saturation'].read_direct(part, selection)
            self.sats[:,0:2] = part / (self.tInt * self.srate)
            part = numpy.empty((self.iDuration, 2), dtype=tuning2['Saturation'].dtype)
            tuning2['Saturation'].read_direct(part, selection)
            self.sats[:,2:4] = part / (self.tInt * self.srate)
            del part
        except KeyError:
            pass
        # Mask out negative saturation values since that indicates the data is
        # not available
        self.sats = numpy.ma.array(self.sats, mask=(self.sats < 0))

        mask1 = tuning1.get('Mask', None)
        mask2 = tuning2.get('Mask', None)

        mask = numpy.zeros(self.spec.shape, dtype=bool)

        for p in data_products:
            if mask1 is not None:
                ## Tuning 1
                ind = 4*(1-1) + mapper[p]
                part = numpy.empty((self.iDuration, self.freq1.size), dtype=mask1[p].dtype)
                mask1[p].read_direct(part, selection)
                mask[ind,:,:] = part.astype(numpy.float32)

                del part

            if mask2 is not None:
                ## Tuning 2
                ind = 4*(2-1) + mapper[p]
                part = numpy.empty((self.iDuration, self.freq2.size), dtype=mask2[p].dtype)
                mask2[p].read_direct(part, selection)
                mask[ind,:,:] = part.astype(numpy.float32)

                del part

        self.spec = numpy.ma.array(self.spec, mask=mask)

        # Construct frequency and time master masks to prevent some masked things from getting unmasked
        self.freqMask = nan_mean(self.spec.mask, axis=1)
        self.timeMask = nan_mean(self.spec.mask, axis=2)
        self.freqMask = numpy.where(self.freqMask >= 0.5, True, False)
        self.timeMask = numpy.where(self.timeMask >= 0.5, True, False)

        # Other data to keep around
        self.timesNPZ = numpy.zeros(obs['time'].shape, dtype=obs['time'].dtype)
        obs['time'].read_direct(self.timesNPZ)
        self.timesNPZRestricted = numpy.zeros(self.spec.shape[1], dtype=obs['time'].dtype)
        obs['time'].read_direct(self.timesNPZRestricted, numpy.s_[self.iOffset:self.iOffset+self.iDuration])
        try:
            fmt = obs['time'].attrs['format']
            scl = obs['time'].attrs['scale']
            try:
                fmt = fmt.decode()
                scl = scl.decode()
            except AttributeError:
                pass

            if fmt != 'unix' or scl != 'utc':
                self.timesNPZ = [AstroTime(*t, format=fmt, scale=scl) for t in self.timesNPZ]
                self.timesNPZ = [t.utc.unix for t in self.timesNPZ]
                self.timesNPZ = numpy.array(self.timesNPZ)

                self.timesNPZRestricted = [AstroTime(*t, format=obs['time'].attrs['format'], scale=obs['time'].attrs['scale']) for t in self.timesNPZRestricted]
                self.timesNPZRestricted = [t.utc.unix for t in self.timesNPZRestricted]
                self.timesNPZRestricted = numpy.array(self.timesNPZRestricted)

            else:
                self.timesNPZ = self.timesNPZ["int"] + self.timesNPZ["frac"]

                self.timesNPZRestricted = self.timesNPZRestricted["int"] + self.timesNPZRestricted["frac"]

        except (KeyError, ValueError):
            pass

        # Deal with the potential for aggregated files
        self.tIntActual = self.tInt
        self.tIntOriginal = self.tInt
        self.filenames = None

        # Gather up the target information
        self.data_products = data_products
        self.target = obs.attrs['TargetName']
        self.ra = obs.attrs['RA']
        self.raUnit = obs.attrs['RA_Units']
        self.dec = obs.attrs['Dec']
        self.decUnit = obs.attrs['Dec_Units']
        self.mode = obs.attrs['TrackingMode']
        self.rbw = obs.attrs['RBW']
        self.rbwUnit = obs.attrs['RBW_Units']

        # Close out the file
        h.close()

        # Compute the bandpass fit
        print(" %6.3f s - Computing bandpass fits" % (time.time() - tStart))
        self.computeBandpass()

        # Find the mean spectra
        print(" %6.3f s - Computing mean spectra" % (time.time() - tStart))
        self.mean = nan_mean(self.spec, axis=1)
        self.meanBandpass = nan_mean(self.specBandpass, axis=1)

        # Set default colobars
        print(" %6.3f s - Setting default colorbar ranges" % (time.time() - tStart))
        self.limits = [None,]*self.spec.shape[0]
        self.limitsBandpass = [None,]*self.spec.shape[0]

        toUse = range(self.spec.shape[2]//10, 9*self.spec.shape[2]//10+1)
        for i in range(self.spec.shape[0]):
            self.limits[i] = findLimits(self.spec[i,:,:], usedB=self.usedB)
        for i in range(self.spec.shape[0]):
            self.limitsBandpass[i] = findLimits(self.specBandpass[i,:,toUse], usedB=self.usedB)

        try:
            self.disconnect()
        except:
            pass

        self.frame.edited = False
        self.frame.setSaveButton()

        print(" %6.3f s - Finished preparing data" % (time.time() - tStart))

    def render(self):
        # Clear the old marks
        self.oldMark = None

        # Clear the old figures
        self.frame.figure1a.clf()
        self.frame.figure1b.clf()
        self.frame.figure1c.clf()
        self.frame.figure2.clf()

        self.connect()

    def computeBandpass(self):
        """
        Compute the bandpass fits.
        """

        if self.bandpassType == 'data':
            self.computeBandpassData()
        else:
            self.computeBandpassInstrumental()

    def computeBandpassData(self):
        """
        Compute data-based bandpass fits.
        """

        meanSpec = nan_median(self.spec, axis=1)

        # Come up with an appropriate smoothing window (wd) and order (od)
        ws = int(round(self.spec.shape[2]/10.0))
        ws = min([41, ws])
        if ws % 2 == 0:
            ws += 1
        od = min([9, ws-2])

        bpm2 = []
        for i in range(self.spec.shape[0]):
            bpm = savitzky_golay(meanSpec[i,:], ws, od, deriv=0)
            bpm = numpy.ma.array(bpm, mask=~numpy.isfinite(bpm))

            if bpm.mean() == 0:
                bpm += 1
            bpm2.append( [bpm / bpm.mean(),] )

        # Apply the bandpass correction
        bpm2 = numpy.array(bpm2)
        self.specBandpass = numpy.ma.array(self.spec.data*1.0, mask=self.spec.mask)
        self.specBandpass.data[...] /= bpm2

        return True

    def computeBandpassInstrumental(self):
        """
        Compute instrument-based bandpass fits.
        """

        def getImpedanceMisMatch(freq):
            freq4, imm4 = stations.lwa1.antennas[0].response(dB=False)

            immIntp = interp1d(freq4, imm4, kind='cubic', bounds_error=False)

            imm = immIntp(freq)
            return imm

        def getARXResponse(freq, filter='full', site=stations.lwa1):
            antennas = site.antennas
            f,r = antennas[0].arx.response(filter='split')
            freq2 = f
            respX2 = numpy.zeros_like(r)
            respY2 = numpy.zeros_like(r)
            for i in range(len(antennas)):
                if antennas[i].combined_status != 33:
                    continue
                f,r = antennas[i].arx.response(filter=filter, dB=False)

                if antennas[i].pol == 0:
                    respX2 += r
                else:
                    respY2 += r
            respX2 /= respX2.max()
            respY2 /= respY2.max()

            respXIntp = interp1d(freq2, respX2, kind='cubic', bounds_error=False)
            respYIntp = interp1d(freq2, respY2, kind='cubic', bounds_error=False)

            respX = respXIntp(freq)
            respY = respYIntp(freq)
            respX = numpy.where(numpy.isfinite(respX), respX, 0)
            respY = numpy.where(numpy.isfinite(respY), respY, 0)
            return respX, respY

        def getDRXResponse(freq, filterCode=7):
            srate = FILTER_CODES[filterCode]
            dpf = dp.drx_filter(sample_rate=srate)

            rDRX = dpf(freq-freq.mean())

            return rDRX

        bpm2 = []
        for i in range(self.spec.shape[0]):
            if i // (self.spec.shape[0]//2) == 0:
                freq = self.freq1
            else:
                freq = self.freq2

            BW = freq.max() - freq.min()
            metric = 1e20
            best = 1
            for fc,fb in FILTER_CODES.items():
                if numpy.abs(BW-fb) < metric:
                    metric = numpy.abs(BW-fb)
                    best = fc

            try:
                rIMM = self._bpCache[('IMM', freq[0], freq[-1], freq.size)]
            except KeyError:
                self._bpCache[('IMM', freq[0], freq[-1], freq.size)] = getImpedanceMisMatch(freq)
                rIMM = self._bpCache[('IMM', freq[0], freq[-1], freq.size)]

            try:
                rARXX, rARXY = self._bpCache[('ARX', freq[0], freq[-1], freq.size)]
            except KeyError:
                self._bpCache[('ARX', freq[0], freq[-1], freq.size)] = getARXResponse(freq, filter=self.arxFilter)
                rARXX, rARXY = self._bpCache[('ARX', freq[0], freq[-1], freq.size)]
            if self.linear and i % (self.spec.shape[0]//2) == 0:
                rARX = rARXX
            elif self.linear and i % (self.spec.shape[0]//2) == 3:
                rARX = rARXY
            else:
                rARX = 0.5*rARXX + 0.5*rARXY

            try:
                rDRX = self._bpCache[('DRX', freq[0], freq[-1], freq.size)]
            except KeyError:
                self._bpCache[('DRX', freq[0], freq[-1], freq.size)] = getDRXResponse(freq, filterCode=best)
                rDRX = self._bpCache[('DRX', freq[0], freq[-1], freq.size)]

            bpm = rIMM/rIMM.max() * rARX/rARX.max() * rDRX/rDRX.max()

            if bpm.mean() == 0:
                bpm += 1
            bpm2.append( [bpm / bpm.mean(),] )

        # Apply the bandpass correction
        bpm2 = numpy.array(bpm2)
        self.specBandpass = numpy.ma.array(self.spec.data*1.0, mask=self.spec.mask)
        self.specBandpass.data[...] /= bpm2

        return True

    def draw(self):
        """
        Draw the waterfall diagram and the total power with time.
        """

        if self.index // (self.spec.shape[0]//2) == 0:
            freq = self.freq1
        else:
            freq = self.freq2

        if self.bandpass:
            spec = self.specBandpass
            limits = self.limitsBandpass
        else:
            spec = self.spec
            limits = self.limits

        # Plot 1(a) - Waterfall
        self.frame.figure1a.clf()
        self.ax1a = self.frame.figure1a.gca()
        if self.usedB:
            m = self.ax1a.imshow(to_dB(spec[self.index,:,:]), interpolation='nearest', extent=(freq[0]/1e6, freq[-1]/1e6, self.time[0], self.time[-1]), origin='lower', cmap=self.cmap, norm=self.norm(limits[self.index][0], limits[self.index][1]))
            try:
                cbar = self.frame.figure1a.colorbar(m, use_gridspec=True)
            except:
                if len(self.frame.figure1a.get_axes()) > 1:
                    self.frame.figure1a.delaxes( self.frame.figure1a.get_axes()[-1] )
                cbar = self.frame.figure1a.colorbar(m, ax=self.ax1a)
            cbar.ax.set_ylabel('PSD [arb. dB]')
        else:
            m = self.ax1a.imshow(spec[self.index,:,:], interpolation='nearest', extent=(freq[0]/1e6, freq[-1]/1e6, self.time[0], self.time[-1]), origin='lower', cmap=self.cmap, norm=self.norm(limits[self.index][0], limits[self.index][1]))
            try:
                cbar = self.frame.figure1a.colorbar(m, use_gridspec=True)
            except TypeError:
                if len(self.frame.figure1a.get_axes()) > 1:
                    self.frame.figure1a.delaxes( self.frame.figure1a.get_axes()[-1] )
                cbar = self.frame.figure1a.colorbar(m, ax=self.ax1a)
            cbar.ax.set_ylabel('PSD [arb. lin.]')
        self.ax1a.axis('auto')
        self.ax1a.set_xlim((freq[0]/1e6, freq[-1]/1e6))
        self.ax1a.set_ylim((self.time[0], self.time[-1]))
        self.ax1a.set_xlabel('Frequency [MHz]')
        self.ax1a.set_ylabel('Elapsed Time - %.3f [s]' % (self.iOffset*self.tInt))
        if self.linear:
            tun = self.index // 4 + 1
            ind = self.index % 4
            mapper = {0: 'XX', 1: 'Re(XY)', 2: 'Im(XY)', 3: 'YY'}
            self.ax1a.set_title('Tuning %i, %s' % (tun, mapper[ind]))
        else:
            tun = self.index // 4 + 1
            ind = self.index % 4
            mapper = {0: 'I', 1: 'Q', 2: 'U', 3: 'V'}
            self.ax1a.set_title('Tuning %i, %s' % (tun, mapper[ind]))

        if self.oldMark is not None:
            self.ax1a.plot([-1e20, 1e20], self.oldMark, color='red')

        self.frame.figure1a.tight_layout()
        self.frame.canvas1a.draw()

        # Plot 1(b) - Saturation Fraction
        self.frame.figure1b.clf()
        self.ax1b = self.frame.figure1b.gca()
        self.ax1b.plot(self.sats[:,2*(tun-1)+0], self.time, linestyle='-', color='blue')
        self.ax1b.plot(self.sats[:,2*(tun-1)+1], self.time, linestyle='-', color='green')
        self.ax1b.set_xlim((-0.05, 1.05))
        self.ax1b.set_ylim((self.time[0], self.time[-1]))
        self.ax1b.set_xlabel('Saturation\nFraction')
        self.ax1b.set_ylabel('Elapsed Time - %.3f [s]' % (self.iOffset*self.tInt))
        self.ax1b.xaxis.set_ticks([0.0, 0.25, 0.5, 0.75, 1.0])
        self.ax1b.xaxis.set_ticklabels(['0', '', '0.5', '', '1'])

        if self.oldMark is not None:
            self.ax1b.plot([0, 1], self.oldMark, color='red')
        self.frame.figure1b.tight_layout()
        self.frame.canvas1b.draw()

        # Plot 1(c) - Drift
        self.drift = spec[:,:,spec.shape[2]//8:7*spec.shape[2]//8].mean(axis=2)

        self.frame.figure1c.clf()
        self.ax1c = self.frame.figure1c.gca()
        if self.usedB:
            z = to_dB(self.drift[self.index,:])
            try:
                self.ax1c.scatter(z, self.time, c=z, marker='x', cmap=self.cmap)
            except ValueError:
                self.ax1c.scatter(z, self.time, c=z, marker='x', cmap=self.cmap, vmin=-1, vmax=1)
            self.ax1c.set_xlabel('Inner 75% Mean Power [arb. dB]')
        else:
            z = self.drift[self.index,:]
            try:
                self.ax1c.scatter(z, self.time, c=z, marker='x', cmap=self.cmap)
            except ValueError:
                self.ax1c.scatter(z, self.time, c=z, marker='x', cmap=self.cmap, vmin=-1, vmax=1)
            self.ax1c.set_xlabel('Inner 75% Mean Power [arb. lin.]')
        self.ax1c.set_ylim((self.time[0], self.time[-1]))
        self.ax1c.set_ylabel('Elapsed Time - %.3f [s]' % (self.iOffset*self.tInt))

        if self.oldMark is not None:
            self.ax1c.plot(self.ax1c.get_xlim(), self.oldMark, color='red')

        self.frame.figure1c.tight_layout()
        self.frame.canvas1c.draw()

    def drawSpectrum(self, clickY, fit=None, fitLabel=None):
        """Get the spectrum at a particular point in time."""

        try:
            dataY = int(round(clickY / self.tInt))
        except TypeError:
            return False

        if self.index // (self.spec.shape[0]//2) == 0:
            freq = self.freq1
        else:
            freq = self.freq2

        if self.bandpass:
            spec = self.specBandpass[self.index,dataY,:]
            medianSpec = self.meanBandpass[self.index,:]
            limits = self.limitsBandpass
        else:
            spec = self.spec[self.index,dataY,:]
            medianSpec = self.mean[self.index,:]
            limits = self.limits

        if self.frame.toolbar.mode == 'zoom rect':
            try:
                oldXlim = self.ax2.get_xlim()
                oldYlim = self.ax2.get_ylim()
            except:
                oldXlim = [freq[0]/1e6, freq[-1]/1e6]
                oldYlim = limits[self.index]
        else:
            oldXlim = [freq[0]/1e6, freq[-1]/1e6]
            oldYlim = limits[self.index]

        self.frame.figure2.clf()
        self.ax2 = self.frame.figure2.gca()
        if self.usedB:
            self.ax2.plot(freq/1e6, to_dB(spec), linestyle=' ', marker='o', label='Current', mec='blue', mfc='None')
            self.ax2.plot(freq/1e6, to_dB(medianSpec), label='Mean', alpha=0.5, color='green')
            self.ax2.set_ylabel('PSD [arb. dB]')
        else:
            self.ax2.plot(freq/1e6, spec, linestyle=' ', marker='o', label='Current', mec='blue', mfc='None')
            self.ax2.plot(freq/1e6, medianSpec, label='Mean', alpha=0.5, color='green')
            self.ax2.set_ylabel('PSD [arb. lin.]')

        if fit is not None:
            if fitLabel is None:
                fitLabel = 'Fit'

            if self.usedB:
                self.ax2.plot(freq/1e6, to_dB(fit), linestyle='--', label=fitLabel, color='orange')
            else:
                self.ax2.plot(freq/1e6, fit, linestyle='--', label=fitLabel, color='orange')

        self.ax2.set_xlim(oldXlim)
        self.ax2.set_ylim(oldYlim)
        self.ax2.legend(loc=0)
        self.ax2.set_xlabel('Frequency [MHz]')

        if self.filenames is None:
            if self.bandpass:
                self.ax2.set_title("%s UTC + bandpass" % datetime.utcfromtimestamp(self.timesNPZRestricted[dataY]))
            else:
                self.ax2.set_title("%s UTC" % datetime.utcfromtimestamp(self.timesNPZRestricted[dataY]))
        else:
            if self.bandpass:
                self.ax2.set_title("%s + bandpass" % self.filenames[dataY])
            else:
                self.ax2.set_title(self.filenames[dataY])

        self.frame.figure2.tight_layout()
        self.frame.canvas2.draw()
        self.spectrumClick = clickY

    def makeMark(self, clickY):
        try:
            dataY = int(round(clickY / self.tInt))
        except TypeError:
            return False

        if self.oldMark is not None:
            try:
                self.ax1a.lines[-1].remove()
            except AttributeError:
                try:
                    del self.ax1a.lines[-1]
                except:
                    pass
            try:
                self.ax1b.lines[-1].remove()
            except AttributeError:
                try:
                    del self.ax1b.lines[-1]
                except:
                    pass
            try:
                self.ax1c.lines[-1].remove()
            except AttributeError:
                try:
                    del self.ax1c.lines[-1]
                except:
                    pass

        oldXSizeA = self.ax1a.get_xlim()
        oldYSizeA = self.ax1a.get_ylim()

        oldXSizeB = self.ax1b.get_xlim()
        oldYSizeB = self.ax1b.get_ylim()

        oldXSizeC = self.ax1c.get_xlim()
        oldYSizeC = self.ax1c.get_ylim()

        self.oldMark = [self.time[dataY],]*2
        self.ax1a.plot([-1e20, 1e20], self.oldMark, color='red')
        self.ax1b.plot([-1e20, 1e20], self.oldMark, color='red')
        self.ax1c.plot([-1e20, 1e20], self.oldMark, color='red')

        self.ax1a.set_xlim(oldXSizeA)
        self.ax1a.set_ylim(oldYSizeA)

        self.ax1b.set_xlim(oldXSizeB)
        self.ax1b.set_ylim(oldYSizeB)

        self.ax1c.set_xlim(oldXSizeC)
        self.ax1c.set_ylim(oldYSizeC)

        self.frame.canvas1a.draw()
        self.frame.canvas1b.draw()
        self.frame.canvas1c.draw()

    def suggestMask(self, index):
        """
        Suggest a mask for the current index.
        """

        nAdjust = {'XX': 1, 'YY': 1, 'XY_real': 2, 'XY_imag': 2,
                   'I': 2, 'Q': 2, 'U': 2, 'V': 2}

        if self.index // (self.spec.shape[0]//2) == 0:
            freq = self.freq1
        else:
            freq = self.freq2
        freqP = freq - freq[len(freq)//2]

        bad = numpy.where( (freqP < -self.srate/2*self.bandpassCut) | (freqP > self.srate/2*self.bandpassCut) )[0]
        for b in bad:
            self.spec.mask[index,:,b] = True
            self.specBandpass.mask[index,:,b] = True
            self.freqMask[index,b] = True

        spec = self.spec.data[index,:,:]
        drift = spec.sum(axis=1)
        coeff = robust.polyfit(numpy.arange(drift.size), drift, self.driftOrder)
        fit = numpy.polyval(coeff, numpy.arange(drift.size))
        rDrift = drift / fit

        try:
            mean = robust.mean(rDrift)
            std  = robust.std(rDrift)
        except ValueError:
            mean = rDrift.mean()
            std  = rDrift.std()

        bad = numpy.where( numpy.abs(rDrift - mean) >= self.driftCut*std )[0]
        for b in bad:
            self.spec.mask[index,b,:] = True
            self.specBandpass.mask[index,b,:] = True
            self.timeMask[index,b] = True

        if self.linear:
            mapper = {0: 'XX', 1: 'XY_real', 2: 'XY_imag', 3: 'YY'}
        else:
            mapper = {0: 'I', 1: 'Q', 2: 'U', 3: 'V'}
        N = self.srate//freq.size*self.tIntActual*nAdjust[mapper[index % 4]]
        kurtosis_arr = numpy.zeros((self.kurtosisSec, self.spec.shape[2]))

        secSize = self.spec.shape[1]//self.kurtosisSec
        for k in range(self.kurtosisSec):
            tStart = k*secSize
            tStop  = (k+1)*secSize

            for j in range(self.spec.shape[2]):
                channel = self.spec.data[index,tStart:tStop,j]
                kurtosis_arr[k,j] = spectral_power(channel, N=N)

        kMean = 1.0
        kStd  = skStd(secSize, N)
        # Correction for some averaged data sets
        if self.filenames is not None:
            kurtosis_arr /= kurtosis_arr.mean()

        bad = numpy.where( numpy.abs(kurtosis_arr - kMean) >= self.kurtosisCut*kStd )
        for k,b in zip(bad[0], bad[1]):
            tStart = k*secSize
            tStop  = (k+1)*secSize

            try:
                for j in range(b-2, b+3):
                    self.spec.mask[index,tStart:tStop,j] = True
                    self.specBandpass.mask[index,tStart:tStop,j] = True
                    self.freqMask[index,j] = True
            except IndexError:
                pass

        self.frame.edited = True
        self.frame.setSaveButton()

        return True

    def resetMask(self, index):
        """
        Reset the specified mask.
        """

        self.spec.mask[index,:,:] = False
        self.specBandpass.mask[index,:,:]= False
        self.timeMask[index,:] = False
        self.freqMask[index,:] = False

        self.frame.edited = True
        self.frame.setSaveButton()

        return True

    def connect(self):
        """
        Connect to all the events we need
        """

        self.cidpress1a = self.frame.figure1a.canvas.mpl_connect('button_press_event', self.on_press1a)
        self.cidpress1b = self.frame.figure1b.canvas.mpl_connect('button_press_event', self.on_press1b)
        self.cidpress1c = self.frame.figure1c.canvas.mpl_connect('button_press_event', self.on_press1c)
        self.cidpress2  = self.frame.figure2.canvas.mpl_connect('button_press_event', self.on_press2)
        self.cidkey2    = self.frame.figure2.canvas.mpl_connect('key_press_event', self.on_key2)
        self.cidmotion  = self.frame.figure1a.canvas.mpl_connect('motion_notify_event', self.on_motion)
        self.frame.toolbar.update = self.on_update

    def on_update(self, *args):
        """
        Override the toolbar 'update' operation.
        """

        self.frame.toolbar.set_history_buttons()

    def on_press1a(self, event):
        """
        On button press we will see if the mouse is over us and store some data
        """

        if event.inaxes:
            clickX = event.xdata
            clickY = event.ydata

            if self.index // (self.spec.shape[0]//2) == 0:
                freq = self.freq1
            else:
                freq = self.freq2

            dataY = int(round(clickY / self.tInt))

            if event.button == 1:
                ## Update the current spectrum
                self.drawSpectrum(clickY)
                self.makeMark(clickY)

            elif event.button == 2:
                ## Unmask
                print("Unmasking %s UTC" % datetime.utcfromtimestamp(self.timesNPZRestricted[dataY]))

                self.spec.mask[self.index, dataY, :] = self.freqMask[self.index,:]
                self.specBandpass.mask[self.index, dataY, :] = self.freqMask[self.index,:]
                self.timeMask[self.index, dataY] = False

                self.draw()
                self.drawSpectrum(clickY)

                self.frame.edited = True
                self.frame.setSaveButton()

            elif event.button == 3:
                ## Mask
                print("Masking %s UTC" % datetime.utcfromtimestamp(self.timesNPZRestricted[dataY]))

                self.spec.mask[self.index, dataY, :] = True
                self.specBandpass.mask[self.index, dataY, :] = True
                self.timeMask[self.index, dataY] = True

                self.draw()
                self.drawSpectrum(clickY)

                self.frame.edited = True
                self.frame.setSaveButton()

            else:
                pass

    def on_press1b(self, event):
        """
        On button press we will see if the mouse is over us and store some data
        """

        if event.inaxes:
            clickX = event.xdata
            clickY = event.ydata

            if self.index // (self.spec.shape[0]//2) == 0:
                freq = self.freq1
            else:
                freq = self.freq2

            dataY = int(round(clickY / self.tInt))

            if event.button == 1:
                ## Update the current spectrum
                self.drawSpectrum(clickY)
                self.makeMark(clickY)

            elif event.button == 2:
                ## Unmask
                print("Unmasking %s UTC" % datetime.utcfromtimestamp(self.timesNPZRestricted[dataY]))

                self.spec.mask[self.index, dataY, :] = self.freqMask[self.index,:]
                self.specBandpass.mask[self.index, dataY, :] = self.freqMask[self.index,:]
                self.timeMask[self.index, dataY] = False

                self.draw()
                self.drawSpectrum(clickY)

                self.frame.edited = True
                self.frame.setSaveButton()

            elif event.button == 3:
                ## Mask
                print("Masking %s UTC" % datetime.utcfromtimestamp(self.timesNPZRestricted[dataY]))

                self.spec.mask[self.index, dataY, :] = True
                self.specBandpass.mask[self.index, dataY, :] = True
                self.timeMask[self.index, dataY] = True

                self.draw()
                self.drawSpectrum(clickY)

                self.frame.edited = True
                self.frame.setSaveButton()

            else:
                pass

    def on_press1c(self, event):
        """
        On button press we will see if the mouse is over us and store some data
        """

        if event.inaxes:
            clickX = event.xdata
            clickY = event.ydata

            if self.index // (self.spec.shape[0]//2) == 0:
                freq = self.freq1
            else:
                freq = self.freq2

            scaleX = self.ax1c.get_xlim()
            rangeX = scaleX[1] - scaleX[0]

            scaleY = self.ax1c.get_ylim()
            rangeY = scaleY[1] - scaleY[0]

            dataY = int(round(clickY / self.tInt))
            lower = dataY - 200
            lower = 0 if lower < 0 else lower
            upper = dataY + 200
            upper = self.drift.shape[0]-1 if upper > self.drift.shape[0]-1 else upper

            if self.usedB:
                d =  ((clickX - to_dB(self.drift.data[self.index,lower:upper]))/rangeX)**2
                d += ((clickY - self.time[lower:upper])/rangeY)**2
                d = numpy.sqrt(d)
                best = numpy.where( d == d.min() )[0][0] + lower
                bestD = d[best - lower]

                print("Clicked at %.3f, %.3f => resolved to entry %i at %.3f, %.3f" % (clickX, clickY, best, to_dB(self.drift.data[self.index, best]), self.time[best]))
            else:
                d =  ((clickX - self.drift.data[self.index,lower:upper])/rangeX)**2
                d += ((clickY - self.time[lower:upper])/rangeY)**2
                d = numpy.sqrt(d)
                best = numpy.where( d == d.min() )[0][0] + lower
                bestD = d[best - lower]

                print("Clicked at %.3f, %.3f => resolved to entry %i at %.3f, %.3f" % (clickX, clickY, best, self.drift.data[self.index, best], self.time[best]))

            if event.button == 1:
                ## Update the current spectrum
                self.drawSpectrum(clickY)
                self.makeMark(clickY)

            elif event.button == 2:
                ## Unmask
                print("Unmasking %s UTC" % datetime.utcfromtimestamp(self.timesNPZRestricted[dataY]))

                self.spec.mask[self.index, best, :] = self.freqMask[self.index,:]
                self.specBandpass.mask[self.index, best, :] = self.freqMask[self.index,:]
                self.timeMask[self.index, best] = False

                self.draw()
                self.drawSpectrum(clickY)

                self.frame.edited = True
                self.frame.setSaveButton()

            elif event.button == 3:
                ## Mask
                print("Masking %s UTC" % datetime.utcfromtimestamp(self.timesNPZRestricted[dataY]))

                self.spec.mask[self.index, best, :] = True
                self.specBandpass.mask[self.index, best, :] = True
                self.timeMask[self.index, best] = True

                self.draw()
                self.drawSpectrum(clickY)

                self.frame.edited = True
                self.frame.setSaveButton()

            else:
                pass

    def on_press2(self, event):
        """
        On button press we will see if the mouse is over us and store some data
        """

        if event.inaxes:
            clickX = event.xdata
            clickY = event.ydata

            if self.index // (self.spec.shape[0]//2) == 0:
                freq = self.freq1
            else:
                freq = self.freq2

            dataX = numpy.where(numpy.abs(clickX-freq/1e6) == (numpy.abs(clickX-freq/1e6).min()))[0][0]

            if event.button == 1:
                ## Nothing right now
                pass

            elif event.button == 2:
                ## Unmask
                print("Unmasking %.3f MHz" % (freq[dataX]/1e6,))

                self.spec.mask[self.index, :, dataX] = self.timeMask[self.index,:]
                self.specBandpass.mask[self.index, :, dataX] = self.timeMask[self.index,:]
                self.freqMask[self.index, dataX] = False

                self.draw()
                self.drawSpectrum(self.spectrumClick)

                self.frame.edited = True
                self.frame.setSaveButton()

            elif event.button == 3:
                ## Mask
                print("Masking %.3f MHz" % (freq[dataX]/1e6,))

                self.spec.mask[self.index, :, dataX] = True
                self.specBandpass.mask[self.index, :, dataX] = True
                self.freqMask[self.index, dataX] = True

                self.draw()
                self.drawSpectrum(self.spectrumClick)

                self.frame.edited = True
                self.frame.setSaveButton()

            else:
                pass

    def on_key2(self, event):
        """
        On key press we will see if the mouse is over us and store some data
        """

        if event.inaxes:
            clickX = event.xdata
            clickY = event.ydata

            if self.index // (self.spec.shape[0]//2) == 0:
                freq = self.freq1
            else:
                freq = self.freq2

            dataX = numpy.where(numpy.abs(clickX-freq/1e6) == (numpy.abs(clickX-freq/1e6).min()))[0][0]
            dataY = int(round(self.spectrumClick / self.tInt))

            if self.bandpass:
                spec = self.specBandpass
            else:
                spec = self.spec

            if event.key == 'h':
                ## Help
                print("Spectrum Window Keys:")
                print("  p - print the information about the underlying point")
                print("  s - print statistics about the current frequency")
                print("  w - write the current spectrum to a text file")
                print("  m - mask the frequency under the cursor")
                print("  u - unmask the frequency under the cursor")
                print("  f - pick a boundary for a power law fit")
                print("  c - clear the current power law fit")
                print("  h - print this help message")

            elif event.key == 'p':
                ## Print
                print("Frequency: %.3f MHz" % (freq[dataX]/1e6,))
                if self.usedB:
                    print("Power: %.3f dB" % to_dB(spec.data[self.index, dataY, dataX]))
                else:
                    print("Power: %.3f" % spec.data[self.index, dataY, dataX])
                print("Flagged? %s" % spec.mask[self.index, dataY, dataX])
                print("===")

            elif event.key == 's':
                ## Statistics
                print("Frequency: %.3f MHz" % (freq[dataX]/1e6,))
                print("Masked Power:")
                print("  Mean: %.3f" % spec[self.index, :, dataX].mean())
                print("  Median: %.3f" % numpy.median(spec[self.index, :, dataX]))
                print("  Std. Dev.: %.3f" % spec[self.index, :, dataX].std())
                print("  Skew: %.3f" % skew(spec[self.index, :, dataX]))
                print("  Kurtosis: %.3f" % kurtosis(spec[self.index, :, dataX]))
                print("===")

            elif event.key == 'w':
                ## Write
                outname = "spectrum.txt"
                while os.path.exists(outname):
                    outpath, outname = os.path.split(outname)
                    outname, ext = os.path.splitext(outname)
                    try:
                        outname, N = outname.rsplit('-', 1)
                    except ValueError:
                        N = "0"
                    try:
                        N = int(N, 10)
                    except ValueError:
                        outname = "%s-%s" % (outname, N)
                        N = 0
                    N += 1
                    outname = "%s%s%s-%i%s" % (outpath, os.path.sep if outpath else '', outname, N, ext)
                print("Saving to '%s'" % outname)

                fn = os.path.basename(self.filename)
                if len(fn) > 33:
                    fn = fn[:30]+"..."
                if self.linear:
                    tun = self.index // 4 + 1
                    ind = self.index % 4
                    mapper = {0: 'XX', 1: 'Re(XY)', 2: 'Im(XY)', 3: 'YY'}
                    pol = mapper[ind]
                else:
                    tun = self.index // 4 + 1
                    ind = self.index % 4
                    mapper = {0: 'I', 1: 'Q', 2: 'U', 3: 'V'}
                    pol = mapper[ind]
                dt = datetime.utcfromtimestamp(self.timesNPZRestricted[dataY])

                beam = self.beam
                srate, sunit = bestFreqUnits(self.srate)
                tInt = self.tInt
                isAggregate = False if self.filenames is None else True
                tIntOrg = self.tIntOriginal
                tIntAct = self.tIntActual

                target = self.target
                if self.raUnit.lower().find('h') != -1:
                    ra = ephem.hours(self.ra*numpy.pi/12.0)
                else:
                    ra = ephem.hours(self.ra*numpy.pi/180.0)
                if self.decUnit.lower().find('h') != -1:
                    dec = ephem.degrees(self.dec*numpy.pi/12.0)
                else:
                    dec = ephem.degrees(self.dec*numpy.pi/180.0)
                mode = self.mode

                rbw = self.rbw
                rbwu = self.rbwUnit

                fh = open(outname, 'w')
                fh.write("#############################################\n")
                fh.write("#                                           #\n")
                fh.write("# Input: %-33s  #\n" % fn)
                fh.write("#                                           #\n")
                fh.write("# Beam: %-2i                                  #\n" % beam)
                fh.write("# RA: %-13s                         #\n" % ra)
                fh.write("# Dec: %-13s                        #\n" % dec)
                fh.write("# Observing Mode: %-10s                #\n" % mode)
                fh.write("# Sample Rate: %-11.4f                  #\n" % srate)
                fh.write("# Sample Rate Units: %-10s             #\n" % sunit)
                fh.write("#                                           #\n")
                fh.write("# Tuning: %-1i                                 #\n" % tun)
                fh.write("# Polarization: %-6s                      #\n" % pol)
                fh.write("# Date/Time: %-26s     #\n" % dt)
                fh.write("# Channel Count: %-10i                 #\n" % freq.size)
                fh.write("# Resolution Bandwidth: %-12.6f        #\n" % rbw)
                fh.write("# Resolution Bandwidth Units: %-10s    #\n" % rbwu)
                fh.write("# Integration Time: %-12.6f            #\n" % tInt)
                fh.write("# Integration Time Units: %-10s        #\n" % "s")
                fh.write("# Aggregated File: %-5s                    #\n" % isAggregate)
                fh.write("#                                           #\n")
                fh.write("# Bandpass Applied: %-5s                   #\n" % self.bandpass)
                fh.write("#                                           #\n")
                fh.write("# Columns:                                  #\n")
                fh.write("#  1. Frequency [Hz]                        #\n")
                fh.write("#  2. Power Spectral Desity [arb. - linear] #\n")
                fh.write("#  3. Masked Channel                        #\n")
                fh.write("#                                           #\n")
                fh.write("#############################################\n")
                for i in range(freq.size):
                    fh.write("%13.5f  %13.6f  %13s\n" % (freq[i], spec.data[self.index,dataY,i], spec.mask[self.index,dataY,i]))
                fh.close()

                print("===")

            elif event.key == 'm':
                ## Mask
                print("Masking %.3f MHz" % (freq[dataX]/1e6,))

                self.spec.mask[self.index, :, dataX] = True
                self.specBandpass.mask[self.index, :, dataX] = True
                self.freqMask[self.index, dataX] = True

                self.draw()
                self.drawSpectrum(self.spectrumClick)

                self.frame.edited = True
                self.frame.setSaveButton()

            elif event.key == 'u':
                ## Unmask
                print("Unmasking %.3f MHz" % (freq[dataX]/1e6,))

                self.spec.mask[self.index, :, dataX] = self.timeMask[self.index,:]
                self.specBandpass.mask[self.index, :, dataX] = self.timeMask[self.index,:]
                self.freqMask[self.index, dataX] = False

                self.draw()
                self.drawSpectrum(self.spectrumClick)

                self.frame.edited = True
                self.frame.setSaveButton()

            elif event.key == 'f':
                ## Power law fit
                self._keyPressCache.append( (self.index*1, dataX, dataY) )

                if len(self._keyPressCache) == 2:
                    (i0,f0,t0), (i1,f1,t1) = self._keyPressCache
                    if i0 != i1 or t0 != t1:
                        del self._keyPressCache[0]
                    else:
                        if f1 < f0:
                            temp = f0
                            f0 = f1
                            f1 = temp
                        print("Fitting from %.3f to %.3f MHz" % (freq[f0]/1e6, freq[f1]/1e6))
                        try:
                            x = freq[f0:f1] / freq[f0]
                            y = spec[t0, i0, f0:f1] / spec[t0, i0, f0]

                            coeff = numpy.polyfit(numpy.log10(x), numpy.log10(y), 1)
                            print("Alpha: %.3f" % coeff[0])
                            fit = 10**(numpy.polyval(coeff, numpy.log10(freq / freq[f0]))) * spec[t0, i0, f0]

                            self.draw()
                            self.drawSpectrum(self.spectrumClick, fit=fit, fitLabel='Alpha: %.3f' % coeff[0])
                        except Exception as e:
                            pass

                        self._keyPressCache = []
                elif len(self._keyPressCache) == 1:
                    print("Move the cursor to the other side of the region to fit push 'f'")

            elif event.key == 'c':
                ## Clear power law fit
                fit = None
                self._keyPressCache = []

                self.draw()
                self.drawSpectrum(self.spectrumClick, fit=fit, fitLabel=None)

            else:
                pass

    def on_motion(self, event):
        """
        On mouse motion display the data value under the cursor
        """

        if event.inaxes:
            clickX = event.xdata
            clickY = event.ydata

            if self.index // (self.spec.shape[0]//2) == 0:
                freq = self.freq1
            else:
                freq = self.freq2

            if self.bandpass:
                spec = self.specBandpass
            else:
                spec = self.spec

            dataX = numpy.where(numpy.abs(clickX-freq/1e6) == (numpy.abs(clickX-freq/1e6).min()))[0][0]
            dataY = numpy.where(numpy.abs(clickY-self.time) == (numpy.abs(clickY-self.time).min()))[0][0]
            if not self.spec.mask[self.index, dataY, dataX]:
                if self.usedB:
                    value = to_dB(spec[self.index, dataY, dataX])
                    self.frame.statusbar.config(text="f=%.4f MHz, t=%.4f s, p=%.2f dB" % (clickX, clickY, value))
                else:
                    value = spec[self.index, dataY, dataX]
                    self.frame.statusbar.config(text="f=%.4f MHz, t=%.4f s, p=%.2f" % (clickX, clickY, value))
            else:
                self.frame.statusbar.config(text="")
        else:
            self.frame.statusbar.config(text="")

    def disconnect(self):
        """
        Disconnect all the stored connection ids
        """

        self.frame.figure1a.canvas.mpl_disconnect(self.cidpress1a)
        self.frame.figure1b.canvas.mpl_disconnect(self.cidpress1b)
        self.frame.figure1c.canvas.mpl_disconnect(self.cidpress1c)
        self.frame.figure2.canvas.mpl_disconnect(self.cidpress2)
        self.frame.figure2.canvas.mpl_disconnect(self.cidkey2)
        self.frame.figure1a.canvas.mpl_disconnect(self.cidmotion)


class MainWindow(tk.Tk):
    def __init__(self):
        super().__init__()

        self.dirname = ''
        self.filename = ''
        self.offset = 0.0
        self.duration = -1
        self.data = None
        self.examineFileButton = None
        self.examineWindow = None

        self.edited = False

        self.title("DRX Waterfall Viewer")
        self.geometry('1000x600')

    def render(self):
        self.initUI()
        self.initEvents()

        self.cAdjust = None
        self._loading_frame = None
        self._resize_after_id = None

    def showLoading(self, message="Loading..."):
        """Show a loading overlay with a message"""
        if self._loading_frame is not None:
            self.hideLoading()

        # Create overlay frame that covers the whole window
        self._loading_frame = tk.Frame(self, bg='#f0f0f0')
        self._loading_frame.place(relx=0, rely=0, relwidth=1, relheight=1)

        # Center container
        container = tk.Frame(self._loading_frame, bg='#f0f0f0')
        container.place(relx=0.5, rely=0.5, anchor=tk.CENTER)

        # Loading message
        self._loading_label = tk.Label(
            container,
            text=message,
            font=('TkDefaultFont', 14),
            bg='#f0f0f0',
            fg='#333333'
        )
        self._loading_label.pack(pady=10)

        # Animated dots
        self._loading_dots = 0
        self._animateLoading()

        self.update()

    def _animateLoading(self):
        """Animate the loading indicator"""
        if self._loading_frame is None:
            return

        self._loading_dots = (self._loading_dots + 1) % 4
        dots = '.' * self._loading_dots
        current_text = self._loading_label.cget('text')
        base_text = current_text.rstrip('.')
        self._loading_label.config(text=base_text + dots)
        self.update_idletasks()

        # Schedule next animation frame
        self._loading_after_id = self.after(300, self._animateLoading)

    def hideLoading(self):
        """Hide the loading overlay"""
        if hasattr(self, '_loading_after_id'):
            self.after_cancel(self._loading_after_id)
        if self._loading_frame is not None:
            self._loading_frame.destroy()
            self._loading_frame = None

    def initUI(self):
        # Status bar at the bottom
        self.statusbar = ttk.Label(self, text="", relief=tk.SUNKEN, anchor=tk.W)
        self.statusbar.pack(side=tk.BOTTOM, fill=tk.X)

        # Menu bar
        menubar = tk.Menu(self)

        # File Menu
        self.fileMenu = tk.Menu(menubar, tearoff=0)
        self.fileMenu.add_command(label='Open', command=self.onOpen, accelerator='Ctrl+O')
        self.fileMenu.add_command(label='Save', command=self.onSave, accelerator='Ctrl+S')
        self.fileMenu.add_command(label='Save As', command=self.onSaveAs)
        self.fileMenu.add_separator()
        self.fileMenu.add_command(label='Exit', command=self.onExit, accelerator='Ctrl+Q')
        menubar.add_cascade(label='File', menu=self.fileMenu)

        # Color Menu
        self.colorMenu = tk.Menu(menubar, tearoff=0)
        self.colorMenu.add_command(label='Auto-scale Colorbar', command=self.onAutoscale)
        self.colorMenu.add_command(label='Adjust Contrast', command=self.onAdjust)

        # Color Map submenu
        self.cmapMenu = tk.Menu(self.colorMenu, tearoff=0)
        self.cmap_var = tk.StringVar(value='jet')
        self.cmapMenu.add_radiobutton(label='Paired', variable=self.cmap_var, value='Paired', command=self.onColorMap)
        self.cmapMenu.add_radiobutton(label='Spectral', variable=self.cmap_var, value='Spectral', command=self.onColorMap)
        self.cmapMenu.add_radiobutton(label='Bone', variable=self.cmap_var, value='bone', command=self.onColorMap)
        self.cmapMenu.add_radiobutton(label='Jet', variable=self.cmap_var, value='jet', command=self.onColorMap)
        self.cmapMenu.add_radiobutton(label='Earth', variable=self.cmap_var, value='gist_earth', command=self.onColorMap)
        self.cmapMenu.add_radiobutton(label='Heat', variable=self.cmap_var, value='gist_heat', command=self.onColorMap)
        self.cmapMenu.add_radiobutton(label='NCAR', variable=self.cmap_var, value='gist_ncar', command=self.onColorMap)
        self.cmapMenu.add_radiobutton(label='Rainbow', variable=self.cmap_var, value='gist_rainbow', command=self.onColorMap)
        self.cmapMenu.add_radiobutton(label='Stern', variable=self.cmap_var, value='gist_stern', command=self.onColorMap)
        self.cmapMenu.add_radiobutton(label='Gray', variable=self.cmap_var, value='gist_gray', command=self.onColorMap)
        self.cmapMenu.add_radiobutton(label='Viridis', variable=self.cmap_var, value='viridis', command=self.onColorMap)
        self.cmapMenu.add_radiobutton(label='Magma', variable=self.cmap_var, value='magma', command=self.onColorMap)
        self.cmapMenu.add_separator()
        self.cmap_invert_var = tk.BooleanVar(value=False)
        self.cmapMenu.add_checkbutton(label='Invert', variable=self.cmap_invert_var, command=self.onColorMap)
        self.colorMenu.add_cascade(label='Color Map', menu=self.cmapMenu)

        # Color Stretch submenu
        self.smapMenu = tk.Menu(self.colorMenu, tearoff=0)
        self.stretch_var = tk.StringVar(value='linear')
        self.smapMenu.add_radiobutton(label='Linear', variable=self.stretch_var, value='linear', command=self.onColorStretch)
        self.smapMenu.add_radiobutton(label='Log', variable=self.stretch_var, value='log', command=self.onColorStretch)
        self.smapMenu.add_radiobutton(label='Square Root', variable=self.stretch_var, value='sqrt', command=self.onColorStretch)
        self.smapMenu.add_radiobutton(label='Squared', variable=self.stretch_var, value='sqrd', command=self.onColorStretch)
        self.smapMenu.add_radiobutton(label='ASinh', variable=self.stretch_var, value='asinh', command=self.onColorStretch)
        self.smapMenu.add_radiobutton(label='Sinh', variable=self.stretch_var, value='sinh', command=self.onColorStretch)
        self.smapMenu.add_radiobutton(label='Histogram Equalization', variable=self.stretch_var, value='histeq', command=self.onColorStretch)
        self.colorMenu.add_cascade(label='Color Stretch', menu=self.smapMenu)
        menubar.add_cascade(label='Color', menu=self.colorMenu)

        # Data Menu
        self.dataMenu = tk.Menu(menubar, tearoff=0)
        self.data_product_var = tk.IntVar(value=0)
        self.dataMenu.add_radiobutton(label='Tuning 1, XX', variable=self.data_product_var, value=0, command=self.onTuning1Product1)
        self.dataMenu.add_radiobutton(label='Tuning 1, Re(XY)', variable=self.data_product_var, value=1, command=self.onTuning1Product2)
        self.dataMenu.add_radiobutton(label='Tuning 1, Im(XY)', variable=self.data_product_var, value=2, command=self.onTuning1Product3)
        self.dataMenu.add_radiobutton(label='Tuning 1, YY', variable=self.data_product_var, value=3, command=self.onTuning1Product4)
        self.dataMenu.add_separator()
        self.dataMenu.add_radiobutton(label='Tuning 2, XX', variable=self.data_product_var, value=4, command=self.onTuning2Product1)
        self.dataMenu.add_radiobutton(label='Tuning 2, Re(XY)', variable=self.data_product_var, value=5, command=self.onTuning2Product2)
        self.dataMenu.add_radiobutton(label='Tuning 2, Im(XY)', variable=self.data_product_var, value=6, command=self.onTuning2Product3)
        self.dataMenu.add_radiobutton(label='Tuning 2, YY', variable=self.data_product_var, value=7, command=self.onTuning2Product4)
        self.dataMenu.add_separator()
        self.dataMenu.add_command(label='Change Time Range', command=self.onRangeChange)
        self.dataMenu.add_command(label='Change Observation', command=self.onObservationChange)
        menubar.add_cascade(label='Data', menu=self.dataMenu)

        # Mask Menu
        self.maskMenu = tk.Menu(menubar, tearoff=0)
        self.maskMenu.add_command(label='Suggest Mask - Current', command=self.onMaskSuggestCurrent)
        self.maskMenu.add_command(label='Suggest Mask - All', command=self.onMaskSuggestAll)
        self.maskMenu.add_separator()
        self.maskMenu.add_command(label='Reset Mask - Current', command=self.onMaskResetCurrent)
        self.maskMenu.add_command(label='Reset Mask - All', command=self.onMaskResetAll)
        self.maskMenu.add_separator()
        self.maskMenu.add_command(label='Adjust Masking Parameters', command=self.onMaskTweak)
        menubar.add_cascade(label='Mask', menu=self.maskMenu)

        # Bandpass Menu
        self.bandpassMenu = tk.Menu(menubar, tearoff=0)
        self.bandpass_var = tk.StringVar(value='off')
        self.bandpassMenu.add_radiobutton(label='Off', variable=self.bandpass_var, value='off', command=self.onBandpassOff)
        self.bandpassMenu.add_radiobutton(label='On', variable=self.bandpass_var, value='on', command=self.onBandpassOn)
        self.bandpassMenu.add_separator()
        self.bandpassMenu.add_command(label='Recompute Fits', command=self.onBandpassRecompute)
        menubar.add_cascade(label='Bandpass', menu=self.bandpassMenu)

        # Details Menu
        self.detailsMenu = tk.Menu(menubar, tearoff=0)
        self.detailsMenu.add_command(label='Current File Info.', command=self.onFileDetails)
        self.detailsMenu.add_command(label='Examine File', command=self.onExamineFile)
        self.detailsMenu.add_separator()
        self.detailsMenu.add_command(label='Zoomable Waterfall', command=self.onZoomableWaterfall)
        self.detailsMenu.add_command(label='Zoomable Drift Curve', command=self.onZoomableDrift)
        self.detailsMenu.add_command(label='Zoomable Power Spectrum', command=self.onZoomablePowerSpectrum)
        menubar.add_cascade(label='Details', menu=self.detailsMenu)

        # Help Menu
        self.helpMenu = tk.Menu(menubar, tearoff=0)
        self.helpMenu.add_command(label='plotHDF Handbook', command=self.onHelp, accelerator='F1')
        self.helpMenu.add_separator()
        self.helpMenu.add_command(label='About', command=self.onAbout)
        menubar.add_cascade(label='Help', menu=self.helpMenu)

        self.config(menu=menubar)

        # Main frame for plots
        main_frame = ttk.Frame(self)
        main_frame.pack(fill=tk.BOTH, expand=True)

        # Top panel with waterfall, saturation, and drift plots
        panel1 = ttk.Frame(main_frame)
        panel1.pack(fill=tk.BOTH, expand=True)

        # Waterfall plot
        self.figure1a = Figure(figsize=(2, 2))
        self.canvas1a = FigureCanvasTkAgg(self.figure1a, master=panel1)
        self.canvas1a.get_tk_widget().pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

        # Saturation fraction plot
        self.figure1b = Figure(figsize=(1, 2))
        self.canvas1b = FigureCanvasTkAgg(self.figure1b, master=panel1)
        self.canvas1b.get_tk_widget().pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

        # Total power with time plot
        self.figure1c = Figure(figsize=(2, 2))
        self.canvas1c = FigureCanvasTkAgg(self.figure1c, master=panel1)
        self.canvas1c.get_tk_widget().pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

        self.panel1 = panel1

        # Bottom panel with spectrum plot and toolbar
        panel3 = ttk.Frame(main_frame)
        panel3.pack(fill=tk.BOTH, expand=True)

        self.figure2 = Figure(figsize=(5, 2))
        self.canvas2 = FigureCanvasTkAgg(self.figure2, master=panel3)
        self.canvas2.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        self.toolbar = NavigationToolbar2Tk(self.canvas2, panel3)
        self.toolbar.update()
        self.toolbar.set_message = lambda s: None  # Disable default coordinate display

        self.panel3 = panel3

    def initEvents(self):
        # Keyboard shortcuts
        self.bind('<Control-o>', lambda e: self.onOpen())
        self.bind('<Control-s>', lambda e: self.onSave())
        self.bind('<Control-q>', lambda e: self.onExit())
        self.bind('<F1>', lambda e: self.onHelp())

        # Arrow key navigation for spectrum
        self.bind('<Up>', lambda e: self.onKeyNavForward())
        self.bind('<Right>', lambda e: self.onKeyNavForward())
        self.bind('<Down>', lambda e: self.onKeyNavBackward())
        self.bind('<Left>', lambda e: self.onKeyNavBackward())

        # Mask/unmask keys
        self.bind('m', lambda e: self.onKeyMask())
        self.bind('M', lambda e: self.onKeyMask())
        self.bind('u', lambda e: self.onKeyUnmask())
        self.bind('U', lambda e: self.onKeyUnmask())

        # Window close
        self.protocol("WM_DELETE_WINDOW", self.onExit)

        # Resize handling
        self.bind('<Configure>', self.onSize)

    def onSize(self, event):
        """Handle window resize with debouncing"""
        # Cancel any pending resize callback to avoid redundant redraws
        if self._resize_after_id is not None:
            self.after_cancel(self._resize_after_id)
        self._resize_after_id = self.after(150, self.resizePlots)

    def resizePlots(self):
        """Resize plots after window resize"""
        self._resize_after_id = None
        try:
            self.figure1a.tight_layout()
            self.canvas1a.draw()
            self.figure1b.tight_layout()
            self.canvas1b.draw()
            self.figure1c.tight_layout()
            self.canvas1c.draw()
            self.figure2.tight_layout()
            self.canvas2.draw()
        except:
            pass

    def onKeyNavForward(self):
        """Move the spectrum display forward by one integration time"""
        if self.data is None or self.data.spectrumClick is None:
            return

        newClick = self.data.spectrumClick + self.data.tInt
        if newClick <= self.data.time[-1]:
            self.data.spectrumClick = newClick
            self.data.drawSpectrum(newClick)
            self.data.makeMark(newClick)

    def onKeyNavBackward(self):
        """Move the spectrum display backward by one integration time"""
        if self.data is None or self.data.spectrumClick is None:
            return

        newClick = self.data.spectrumClick - self.data.tInt
        if newClick >= self.data.time[0]:
            self.data.spectrumClick = newClick
            self.data.drawSpectrum(newClick)
            self.data.makeMark(newClick)

    def onKeyMask(self):
        """Mask the current integration"""
        if self.data is None or self.data.spectrumClick is None:
            return

        dataY = int(round(self.data.spectrumClick / self.data.tInt))

        self.data.spec.mask[self.data.index, dataY, :] = True
        self.data.specBandpass.mask[self.data.index, dataY, :] = True
        self.data.timeMask[self.data.index, dataY] = True

        self.data.draw()
        self.data.drawSpectrum(self.data.spectrumClick)

        self.edited = True
        self.setSaveButton()

    def onKeyUnmask(self):
        """Unmask the current integration"""
        if self.data is None or self.data.spectrumClick is None:
            return

        dataY = int(round(self.data.spectrumClick / self.data.tInt))

        self.data.spec.mask[self.data.index, dataY, :] = self.data.freqMask[self.data.index, :]
        self.data.specBandpass.mask[self.data.index, dataY, :] = self.data.freqMask[self.data.index, :]
        self.data.timeMask[self.data.index, dataY] = False

        self.data.draw()
        self.data.drawSpectrum(self.data.spectrumClick)

        self.edited = True
        self.setSaveButton()

    def setSaveButton(self):
        """Update Save menu item state based on edited flag"""
        if self.edited:
            self.fileMenu.entryconfig(1, state='normal')  # Save
        else:
            self.fileMenu.entryconfig(1, state='disabled')  # Save

    def setDataMenuOptions(self):
        """Enable/disable data menu items based on available data products"""
        if self.data is None:
            return

        # First, disable all (skip index 4 which is a separator)
        for i in [0, 1, 2, 3, 5, 6, 7, 8]:
            self.dataMenu.entryconfig(i, state='disabled')

        # Enable based on data products
        mapper = {'XX': 0, 'I': 0, 'XY_real': 1, 'Q': 1, 'XY_imag': 2, 'U': 2, 'YY': 3, 'V': 3}
        for dp in self.data.data_products:
            idx = mapper.get(dp, -1)
            if idx >= 0:
                # Tuning 1
                self.dataMenu.entryconfig(idx, state='normal')
                # Tuning 2
                self.dataMenu.entryconfig(idx + 5, state='normal')  # +5 accounts for separator

        # Set labels based on linear/Stokes
        if self.data.linear:
            labels = ['Tuning 1, XX', 'Tuning 1, Re(XY)', 'Tuning 1, Im(XY)', 'Tuning 1, YY',
                      'Tuning 2, XX', 'Tuning 2, Re(XY)', 'Tuning 2, Im(XY)', 'Tuning 2, YY']
        else:
            labels = ['Tuning 1, I', 'Tuning 1, Q', 'Tuning 1, U', 'Tuning 1, V',
                      'Tuning 2, I', 'Tuning 2, Q', 'Tuning 2, U', 'Tuning 2, V']

        for i, label in enumerate(labels):
            idx = i if i < 4 else i + 1  # Account for separator
            self.dataMenu.entryconfig(idx, label=label)

    # File menu handlers
    def onOpen(self):
        """Open a HDF5 file"""
        filename = filedialog.askopenfilename(
            initialdir=self.dirname,
            title="Select HDF5 Waterfall File",
            filetypes=[("HDF5 files", "*.hdf5"), ("All files", "*.*")]
        )
        if filename:
            self.dirname = os.path.dirname(filename)
            self.filename = filename

            self.config(cursor='watch')
            self.showLoading("Loading " + os.path.basename(filename))

            try:
                self.data.loadData(filename, obsID=1)
                self.data.render()
                self.data.draw()

                # Enable menus
                self._enableMenus(True)
                self.setDataMenuOptions()

                self.edited = False
                self.setSaveButton()
            except Exception as e:
                messagebox.showerror("Error", f"Failed to load file: {e}")

            self.hideLoading()
            self.config(cursor='')

    def onSave(self):
        """Save mask to current file"""
        if not self.filename:
            return

        self.config(cursor='watch')
        self.update()

        try:
            self._saveMask(self.filename)
            self.edited = False
            self.setSaveButton()
        except Exception as e:
            messagebox.showerror("Error", f"Failed to save: {e}")

        self.config(cursor='')

    def onSaveAs(self):
        """Save mask to a new file"""
        filename = filedialog.asksaveasfilename(
            initialdir=self.dirname,
            title="Save HDF5 File As",
            filetypes=[("HDF5 files", "*.hdf5"), ("All files", "*.*")],
            defaultextension=".hdf5"
        )
        if filename:
            self.config(cursor='watch')
            self.update()

            try:
                # Copy file first if different
                if filename != self.filename:
                    import shutil
                    shutil.copy2(self.filename, filename)

                self._saveMask(filename)
                self.filename = filename
                self.dirname = os.path.dirname(filename)
                self.edited = False
                self.setSaveButton()
            except Exception as e:
                messagebox.showerror("Error", f"Failed to save: {e}")

            self.config(cursor='')

    def _saveMask(self, filename):
        """Save mask data to HDF5 file"""
        h = h5py.File(filename, 'a')
        obs = h.get('Observation%i' % self.data.obsID, None)

        tuning1 = obs.get('Tuning1', None)
        tuning2 = obs.get('Tuning2', None)
        if tuning2 is None:
            tuning2 = tuning1

        mask1 = tuning1.get('Mask', None)
        mask2 = tuning2.get('Mask', None)

        mapper = {'XX': 0, 'I': 0, 'XY_real': 1, 'Q': 1, 'XY_imag': 2, 'U': 2, 'YY': 3, 'V': 3}

        for p in self.data.data_products:
            ind1 = 4*(1-1) + mapper[p]
            ind2 = 4*(2-1) + mapper[p]

            selection = numpy.s_[self.data.iOffset:self.data.iOffset+self.data.iDuration, :]

            if mask1 is not None:
                mask1[p].write_direct(self.data.spec.mask[ind1,:,:].astype(mask1[p].dtype), dest_sel=selection)
            if mask2 is not None:
                mask2[p].write_direct(self.data.spec.mask[ind2,:,:].astype(mask2[p].dtype), dest_sel=selection)

        h.close()

    def onExit(self):
        """Exit the application"""
        if self.edited:
            result = messagebox.askyesnocancel("Unsaved Changes",
                "You have unsaved changes. Do you want to save before exiting?")
            if result is None:  # Cancel
                return
            elif result:  # Yes
                self.onSave()

        self.destroy()

    def _enableMenus(self, enable):
        """Enable or disable menus"""
        state = 'normal' if enable else 'disabled'

        # Color menu
        for i in range(self.colorMenu.index('end') + 1):
            try:
                self.colorMenu.entryconfig(i, state=state)
            except:
                pass

        # Data menu
        for i in range(self.dataMenu.index('end') + 1):
            try:
                self.dataMenu.entryconfig(i, state=state)
            except:
                pass

        # Mask menu
        for i in range(self.maskMenu.index('end') + 1):
            try:
                self.maskMenu.entryconfig(i, state=state)
            except:
                pass

        # Bandpass menu
        for i in range(self.bandpassMenu.index('end') + 1):
            try:
                self.bandpassMenu.entryconfig(i, state=state)
            except:
                pass

        # Details menu
        for i in range(self.detailsMenu.index('end') + 1):
            try:
                self.detailsMenu.entryconfig(i, state=state)
            except:
                pass

    # Color menu handlers
    def onAutoscale(self):
        """Auto-scale the colorbar using 5th/99th percentile"""
        if self.data is None:
            return

        self.config(cursor='watch')
        self.update()

        i = self.data.index
        toUse = numpy.arange(self.data.spec.shape[2]//10, 9*self.data.spec.shape[2]//10)

        if self.data.bandpass:
            self.data.limitsBandpass[i] = list(nan_percentile(self.data.specBandpass[i,:,toUse], (5, 99)))
        else:
            self.data.limits[i] = list(nan_percentile(self.data.spec[i,:,:], (5, 99)))

        # Convert to dB if needed
        if self.data.usedB:
            if self.data.bandpass:
                self.data.limitsBandpass[i] = [to_dB(v) for v in self.data.limitsBandpass[i]]
            else:
                self.data.limits[i] = [to_dB(v) for v in self.data.limits[i]]

        # Handle NaN/infinity cases
        if self.data.bandpass:
            if not numpy.isfinite(self.data.limitsBandpass[i][0]) and not numpy.isfinite(self.data.limitsBandpass[i][1]):
                self.data.limitsBandpass[i] = [0, 1]
            elif not numpy.isfinite(self.data.limitsBandpass[i][0]):
                self.data.limitsBandpass[i] = [self.data.limitsBandpass[i][1]+1, self.data.limitsBandpass[i][1]]
            elif not numpy.isfinite(self.data.limitsBandpass[i][1]):
                self.data.limitsBandpass[i] = [self.data.limitsBandpass[i][0], self.data.limitsBandpass[i][0]-1]
        else:
            if not numpy.isfinite(self.data.limits[i][0]) and not numpy.isfinite(self.data.limits[i][1]):
                self.data.limits[i] = [0, 1]
            elif not numpy.isfinite(self.data.limits[i][0]):
                self.data.limits[i] = [self.data.limits[i][1]+1, self.data.limits[i][1]]
            elif not numpy.isfinite(self.data.limits[i][1]):
                self.data.limits[i] = [self.data.limits[i][0], self.data.limits[i][0]-1]

        self.data.draw()
        self.data.drawSpectrum(self.data.spectrumClick)
        self.data.makeMark(self.data.spectrumClick)

        self.config(cursor='')

    def onAdjust(self):
        """Open contrast adjustment dialog"""
        if self.data is None:
            return

        if self.cAdjust is not None:
            self.cAdjust.lift()
            return

        self.cAdjust = ContrastAdjust(self)

    def onColorMap(self):
        """Change the color map"""
        if self.data is None:
            return

        cmap_name = self.cmap_var.get()
        if self.cmap_invert_var.get():
            cmap_name = cmap_name + '_r'

        self.data.cmap = cm.get_cmap(cmap_name)
        self.data.draw()

    def onColorStretch(self):
        """Change the color stretch"""
        if self.data is None:
            return

        stretch = self.stretch_var.get()
        stretch_map = {
            'linear': Normalize,
            'log': LogNorm,
            'sqrt': SqrtNorm,
            'sqrd': SqrdNorm,
            'asinh': AsinhNorm,
            'sinh': SinhNorm,
            'histeq': HistEqNorm
        }
        self.data.norm = stretch_map.get(stretch, Normalize)
        self.data.draw()

    # Data menu handlers
    def _switchProduct(self, index):
        """Switch to a different data product"""
        if self.data is None:
            return

        self.config(cursor='watch')
        self.update()

        self.data.index = index
        self.data.draw()
        if self.data.spectrumClick is not None:
            self.data.drawSpectrum(self.data.spectrumClick)

        self.config(cursor='')

    def onTuning1Product1(self):
        self._switchProduct(0)

    def onTuning1Product2(self):
        self._switchProduct(1)

    def onTuning1Product3(self):
        self._switchProduct(2)

    def onTuning1Product4(self):
        self._switchProduct(3)

    def onTuning2Product1(self):
        self._switchProduct(4)

    def onTuning2Product2(self):
        self._switchProduct(5)

    def onTuning2Product3(self):
        self._switchProduct(6)

    def onTuning2Product4(self):
        self._switchProduct(7)

    def onRangeChange(self):
        """Open time range adjustment dialog"""
        if self.data is None:
            return
        TimeRangeAdjust(self, mode='Adjust')

    def onObservationChange(self):
        """Open observation change dialog"""
        if self.data is None:
            return
        SwitchObservation(self)

    # Mask menu handlers
    def onMaskSuggestCurrent(self):
        """Suggest mask for current data product"""
        if self.data is None:
            return

        self.config(cursor='watch')
        self.update()

        self.data.suggestMask(self.data.index)
        self.data.draw()
        self.data.drawSpectrum(self.data.spectrumClick)
        self.data.makeMark(self.data.spectrumClick)

        self.config(cursor='')

    def onMaskSuggestAll(self):
        """Suggest mask for all data products"""
        if self.data is None:
            return

        self.config(cursor='watch')
        self.update()

        for i in range(self.data.spec.shape[0]):
            self.data.suggestMask(i)
        self.data.draw()
        self.data.drawSpectrum(self.data.spectrumClick)
        self.data.makeMark(self.data.spectrumClick)

        self.config(cursor='')

    def onMaskResetCurrent(self):
        """Reset mask for current data product"""
        if self.data is None:
            return

        self.config(cursor='watch')
        self.update()

        self.data.resetMask(self.data.index)
        self.data.draw()
        self.data.drawSpectrum(self.data.spectrumClick)
        self.data.makeMark(self.data.spectrumClick)

        self.config(cursor='')

    def onMaskResetAll(self):
        """Reset mask for all data products"""
        if self.data is None:
            return

        self.config(cursor='watch')
        self.update()

        for i in range(self.data.spec.shape[0]):
            self.data.resetMask(i)
        self.data.draw()
        self.data.drawSpectrum(self.data.spectrumClick)
        self.data.makeMark(self.data.spectrumClick)

        self.config(cursor='')

    def onMaskTweak(self):
        """Open masking parameters dialog"""
        if self.data is None:
            return
        MaskingAdjust(self)

    # Bandpass menu handlers
    def onBandpassOn(self):
        """Turn on bandpass correction"""
        if self.data is None:
            return

        self.data.bandpass = True
        self.data.draw()
        self.data.drawSpectrum(self.data.spectrumClick)
        self.data.makeMark(self.data.spectrumClick)

    def onBandpassOff(self):
        """Turn off bandpass correction"""
        if self.data is None:
            return

        self.data.bandpass = False
        self.data.draw()
        self.data.drawSpectrum(self.data.spectrumClick)
        self.data.makeMark(self.data.spectrumClick)

    def onBandpassRecompute(self):
        """Recompute bandpass fits"""
        if self.data is None:
            return

        self.config(cursor='watch')
        self.update()

        self.data.computeBandpass()
        self.data.meanBandpass = nan_mean(self.data.specBandpass, axis=1)

        toUse = range(self.data.spec.shape[2]//10, 9*self.data.spec.shape[2]//10+1)
        for i in range(self.data.spec.shape[0]):
            self.data.limitsBandpass[i] = findLimits(self.data.specBandpass[i,:,toUse], usedB=self.data.usedB)

        self.data.draw()
        if self.data.bandpass:
            self.data.drawSpectrum(self.data.spectrumClick)
            self.data.makeMark(self.data.spectrumClick)
        self.config(cursor='')

    # Details menu handlers
    def onFileDetails(self):
        """Show file details"""
        if self.data is None:
            return

        beam = self.data.beam
        srate, sunit = bestFreqUnits(self.data.srate)
        tInt = self.data.tInt

        target = self.data.target
        if self.data.raUnit.lower().find('h') != -1:
            ra = ephem.hours(self.data.ra*numpy.pi/12.0)
        else:
            ra = ephem.hours(self.data.ra*numpy.pi/180.0)
        if self.data.decUnit.lower().find('h') != -1:
            dec = ephem.degrees(self.data.dec*numpy.pi/12.0)
        else:
            dec = ephem.degrees(self.data.dec*numpy.pi/180.0)
        mode = self.data.mode

        rbw = self.data.rbw
        rbwu = self.data.rbwUnit

        info = f"""File: {os.path.basename(self.filename)}

Target: {target}
RA: {ra}
Dec: {dec}
Tracking Mode: {mode}

Beam: {beam}
Sample Rate: {srate:.4f} {sunit}
Integration Time: {tInt:.6f} s
Resolution Bandwidth: {rbw:.6f} {rbwu}

Channels: {self.data.freq1.size}
Integrations: {self.data.spec.shape[1]}
"""

        messagebox.showinfo("File Information", info)

    def onExamineFile(self):
        """
        For aggregated data sets, open a new plotHDF window to examine the
        sub-file at the current time selection.
        """
        if self.data is None or self.data.filenames is None:
            messagebox.showinfo("Info", "No sub-files to examine")
            return

        # Check if another examine window is already running
        if self.examineWindow is not None:
            self.examineWindow.poll()
            if self.examineWindow.returncode is None:
                messagebox.showwarning("Warning",
                    "Another sub-file examination window is already running")
                return
            else:
                self.examineWindow = None

        # Make sure we have clicked to select a time
        if self.data.spectrumClick is None:
            messagebox.showinfo("Info", "Please click on the waterfall to select a time first")
            return

        try:
            dataY = int(round(self.data.spectrumClick / self.data.tInt))
        except:
            messagebox.showinfo("Info", "No sub-file currently selected")
            return

        # Make sure dataY is within bounds
        if dataY < 0 or dataY >= len(self.data.filenames):
            messagebox.showerror("Error", "Selected time is out of range")
            return

        # Get the filename and check if it exists
        filename = self.data.filenames[dataY]
        if not os.path.exists(filename):
            # Try the path relative to the current file
            basepath = os.path.dirname(self.data.filename)
            basename = os.path.basename(self.data.filenames[dataY])
            filename = os.path.join(basepath, basename)

            if not os.path.exists(filename):
                messagebox.showerror("Error", f"Cannot find file: {filename}")
                return

        # Get the current script path
        script = sys.argv[0]

        # Launch new instance
        self.examineWindow = subprocess.Popen([sys.executable, script, filename])

    def onZoomableWaterfall(self):
        """Open zoomable waterfall window"""
        if self.data is None:
            return
        WaterfallDisplay(self)

    def onZoomableDrift(self):
        """Open zoomable drift curve window"""
        if self.data is None:
            return
        DriftCurveDisplay(self)

    def onZoomablePowerSpectrum(self):
        """Open zoomable power spectrum window"""
        if self.data is None:
            return
        PowerSpectrumDisplay(self)

    # Help menu handlers
    def onHelp(self):
        """Show help window"""
        HelpWindow(self)

    def onAbout(self):
        """Show about dialog"""
        about_text = f"""plotHDF - DRX HDF5 Waterfall Viewer

Version: {__version__}
Author: {__author__}

A tool for viewing and editing DRX waterfall HDF5 files.

Tkinter version converted from wxPython."""

        messagebox.showinfo("About plotHDF", about_text)

    def getToolBar(self):
        return self.toolbar


class ContrastAdjust(tk.Toplevel):
    """Dialog for adjusting contrast limits"""

    def __init__(self, parent):
        super().__init__(parent)
        self.parent = parent

        self.title("Contrast Adjustment")
        self.transient(parent)

        self.protocol("WM_DELETE_WINDOW", self.onClose)

        self._build_ui()
        self.parent.cAdjust = self

        # Bind Enter key to apply
        self.bind('<Return>', self.onKeyPress)

    def _build_ui(self):
        frame = ttk.Frame(self, padding=10)
        frame.pack(fill=tk.BOTH, expand=True)

        # Get current limits
        data = self.parent.data
        if data.bandpass:
            limits = data.limitsBandpass[data.index]
            label_text = 'Tuning %i, Pol. %i - Bandpass' % (data.index//2+1, data.index%2)
        else:
            limits = data.limits[data.index]
            label_text = 'Tuning %i, Pol. %i' % (data.index//2+1, data.index%2)

        # Title label
        ttk.Label(frame, text=label_text).grid(row=0, column=0, columnspan=4, sticky=tk.W, pady=5)

        # Upper Limit
        ttk.Label(frame, text="Upper Limit:").grid(row=1, column=0, sticky=tk.W, pady=5)
        self.upper_var = tk.StringVar(value='%.1f' % limits[1])
        self.upper_entry = ttk.Entry(frame, textvariable=self.upper_var, width=12)
        self.upper_entry.grid(row=1, column=1, pady=5, padx=5)
        ttk.Button(frame, text="-", width=3, command=self.onUpperDecrease).grid(row=1, column=2, pady=5, padx=2)
        ttk.Button(frame, text="+", width=3, command=self.onUpperIncrease).grid(row=1, column=3, pady=5, padx=2)

        # Lower Limit
        ttk.Label(frame, text="Lower Limit:").grid(row=2, column=0, sticky=tk.W, pady=5)
        self.lower_var = tk.StringVar(value='%.1f' % limits[0])
        self.lower_entry = ttk.Entry(frame, textvariable=self.lower_var, width=12)
        self.lower_entry.grid(row=2, column=1, pady=5, padx=5)
        ttk.Button(frame, text="-", width=3, command=self.onLowerDecrease).grid(row=2, column=2, pady=5, padx=2)
        ttk.Button(frame, text="+", width=3, command=self.onLowerIncrease).grid(row=2, column=3, pady=5, padx=2)

        # Range (read-only)
        ttk.Label(frame, text="Range:").grid(row=3, column=0, sticky=tk.W, pady=5)
        self.range_var = tk.StringVar(value='%.1f' % self._getRange())
        self.range_entry = ttk.Entry(frame, textvariable=self.range_var, width=12, state='readonly')
        self.range_entry.grid(row=3, column=1, pady=5, padx=5)

        # Separator
        ttk.Separator(frame, orient=tk.HORIZONTAL).grid(row=4, column=0, columnspan=4, sticky='ew', pady=10)

        # Ok button
        ttk.Button(frame, text="Ok", command=self.onOk).grid(row=5, column=3, pady=5)

    def _getRange(self):
        data = self.parent.data
        index = data.index
        if data.bandpass:
            return data.limitsBandpass[index][1] - data.limitsBandpass[index][0]
        else:
            return data.limits[index][1] - data.limits[index][0]

    def _getIncrement(self):
        return 0.1 * self._getRange()

    def _updateRange(self):
        self.range_var.set('%.1f' % self._getRange())

    def onUpperDecrease(self):
        data = self.parent.data
        index = data.index
        increment = self._getIncrement()
        if data.bandpass:
            data.limitsBandpass[index][1] -= increment
            self.upper_var.set('%.1f' % data.limitsBandpass[index][1])
        else:
            data.limits[index][1] -= increment
            self.upper_var.set('%.1f' % data.limits[index][1])
        self._updateRange()
        data.draw()

    def onUpperIncrease(self):
        data = self.parent.data
        index = data.index
        increment = self._getIncrement()
        if data.bandpass:
            data.limitsBandpass[index][1] += increment
            self.upper_var.set('%.1f' % data.limitsBandpass[index][1])
        else:
            data.limits[index][1] += increment
            self.upper_var.set('%.1f' % data.limits[index][1])
        self._updateRange()
        data.draw()

    def onLowerDecrease(self):
        data = self.parent.data
        index = data.index
        increment = self._getIncrement()
        if data.bandpass:
            data.limitsBandpass[index][0] -= increment
            self.lower_var.set('%.1f' % data.limitsBandpass[index][0])
        else:
            data.limits[index][0] -= increment
            self.lower_var.set('%.1f' % data.limits[index][0])
        self._updateRange()
        data.draw()

    def onLowerIncrease(self):
        data = self.parent.data
        index = data.index
        increment = self._getIncrement()
        if data.bandpass:
            data.limitsBandpass[index][0] += increment
            self.lower_var.set('%.1f' % data.limitsBandpass[index][0])
        else:
            data.limits[index][0] += increment
            self.lower_var.set('%.1f' % data.limits[index][0])
        self._updateRange()
        data.draw()

    def onKeyPress(self, event):
        """Apply values when Enter is pressed"""
        try:
            data = self.parent.data
            index = data.index
            lower_val = float(self.lower_var.get())
            upper_val = float(self.upper_var.get())
            if data.bandpass:
                data.limitsBandpass[index][0] = lower_val
                data.limitsBandpass[index][1] = upper_val
            else:
                data.limits[index][0] = lower_val
                data.limits[index][1] = upper_val
            self._updateRange()
            data.draw()
        except ValueError:
            pass

    def onOk(self):
        try:
            data = self.parent.data
            index = data.index
            lower_val = float(self.lower_var.get())
            upper_val = float(self.upper_var.get())
            if data.bandpass:
                data.limitsBandpass[index][0] = lower_val
                data.limitsBandpass[index][1] = upper_val
            else:
                data.limits[index][0] = lower_val
                data.limits[index][1] = upper_val
            data.draw()
        except ValueError:
            messagebox.showerror("Error", "Please enter valid numbers")
            return

        self.parent.cAdjust = None
        self.destroy()

    def onClose(self):
        self.parent.cAdjust = None
        self.destroy()


class TimeRangeAdjust(tk.Toplevel):
    """Dialog for adjusting time range"""

    def __init__(self, parent, mode='Adjust'):
        super().__init__(parent)
        self.parent = parent
        self.mode = mode

        self.title("Adjust Time Range")
        self.transient(parent)
        self.geometry('300x180')

        self._build_ui()

    def _build_ui(self):
        frame = ttk.Frame(self, padding=10)
        frame.pack(fill=tk.BOTH, expand=True)

        # Offset
        ttk.Label(frame, text="Offset (s):").grid(row=0, column=0, sticky=tk.W, pady=5)
        self.offset_var = tk.StringVar(value=str(self.parent.offset))
        self.offset_entry = ttk.Entry(frame, textvariable=self.offset_var, width=15)
        self.offset_entry.grid(row=0, column=1, pady=5)

        # Duration
        ttk.Label(frame, text="Duration (s):").grid(row=1, column=0, sticky=tk.W, pady=5)
        dur = self.parent.duration if self.parent.duration > 0 else -1
        self.duration_var = tk.StringVar(value=str(dur))
        self.duration_entry = ttk.Entry(frame, textvariable=self.duration_var, width=15)
        self.duration_entry.grid(row=1, column=1, pady=5)

        ttk.Label(frame, text="(-1 for all)").grid(row=1, column=2, sticky=tk.W, pady=5)

        # Whole observation checkbox
        self.whole_obs_var = tk.BooleanVar(value=(self.parent.offset == 0 and self.parent.duration < 0))
        self.whole_obs_check = ttk.Checkbutton(
            frame, text="Display whole observation",
            variable=self.whole_obs_var, command=self.onWholeObsToggle
        )
        self.whole_obs_check.grid(row=2, column=0, columnspan=3, sticky=tk.W, pady=5)

        # Set initial state of entries based on checkbox
        if self.whole_obs_var.get():
            self.offset_entry.config(state='disabled')
            self.duration_entry.config(state='disabled')

        # Buttons
        btn_frame = ttk.Frame(frame)
        btn_frame.grid(row=3, column=0, columnspan=3, pady=10)

        ttk.Button(btn_frame, text="Apply", command=self.onApply).pack(side=tk.LEFT, padx=5)
        cancel_btn = ttk.Button(btn_frame, text="Cancel", command=self.destroy)
        cancel_btn.pack(side=tk.LEFT, padx=5)

        # Disable cancel in 'New' mode
        if self.mode != 'Adjust':
            cancel_btn.config(state='disabled')

    def onWholeObsToggle(self):
        """Enable/disable offset and duration entries based on checkbox"""
        if self.whole_obs_var.get():
            self.offset_entry.config(state='disabled')
            self.duration_entry.config(state='disabled')
        else:
            self.offset_entry.config(state='normal')
            self.duration_entry.config(state='normal')

    def onApply(self):
        try:
            if self.whole_obs_var.get():
                offset = 0.0
                duration = -1.0
            else:
                offset = float(self.offset_var.get())
                duration = float(self.duration_var.get())

            self.parent.offset = offset
            self.parent.duration = duration

            # Reload data
            self.parent.config(cursor='watch')
            self.parent.update()

            self.parent.data.loadData(self.parent.filename, obsID=self.parent.data.obsID)
            self.parent.data.render()
            self.parent.data.draw()

            self.parent.config(cursor='')
            self.destroy()
        except ValueError:
            messagebox.showerror("Error", "Please enter valid numbers")


class SwitchObservation(tk.Toplevel):
    """Dialog for switching observations"""

    def __init__(self, parent):
        super().__init__(parent)
        self.parent = parent

        self.title("Switch Observation")
        self.transient(parent)

        self._build_ui()

    def _build_ui(self):
        frame = ttk.Frame(self, padding=10)
        frame.pack(fill=tk.BOTH, expand=True)

        # Get observation list with details
        self.obs_list = self._getObsList()

        # Radio button variable
        self.obs_var = tk.IntVar(value=self.parent.data.obsID)

        # Create radio buttons for each observation
        for i, (obs_id, target_name, tracking_mode) in enumerate(self.obs_list):
            label = f"{target_name} in mode {tracking_mode} (#{obs_id})"
            rb = ttk.Radiobutton(frame, text=label, variable=self.obs_var, value=obs_id)
            rb.grid(row=i, column=0, sticky=tk.W, pady=2)

        # Buttons
        btn_frame = ttk.Frame(frame)
        btn_frame.grid(row=len(self.obs_list), column=0, pady=10)

        ttk.Button(btn_frame, text="Ok", command=self.onSwitch).pack(side=tk.LEFT, padx=5)
        ttk.Button(btn_frame, text="Cancel", command=self.destroy).pack(side=tk.LEFT, padx=5)

    def _getObsList(self):
        """Return a list of tuples (obsID, TargetName, TrackingMode) for each observation"""
        obs_list = []
        h = h5py.File(self.parent.data.filename, 'r')
        for key in sorted(h.keys()):
            if key.startswith('Observation'):
                obs_id = int(key[11:])
                obs = h.get(key, None)
                target_name = obs.attrs.get('TargetName', 'Unknown')
                if isinstance(target_name, bytes):
                    target_name = target_name.decode('utf-8')
                tracking_mode = obs.attrs.get('TrackingMode', 'Unknown')
                if isinstance(tracking_mode, bytes):
                    tracking_mode = tracking_mode.decode('utf-8')
                obs_list.append((obs_id, target_name, tracking_mode))
        h.close()
        return obs_list

    def onSwitch(self):
        obs_id = self.obs_var.get()

        if obs_id == self.parent.data.obsID:
            # No change needed
            self.destroy()
            return

        self.parent.config(cursor='watch')
        self.parent.update()

        try:
            self.parent.data.loadData(self.parent.filename, obsID=obs_id)
            self.parent.data.render()
            self.parent.data.draw()
        except Exception as e:
            messagebox.showerror("Error", f"Failed to switch observation: {e}")

        self.parent.config(cursor='')
        self.destroy()


class MaskingAdjust(tk.Toplevel):
    """Dialog for adjusting masking parameters"""

    def __init__(self, parent):
        super().__init__(parent)
        self.parent = parent

        self.title("Masking Parameters")
        self.transient(parent)

        self._build_ui()

        self.parent.cAdjust = self

    def _build_ui(self):
        frame = ttk.Frame(self, padding=10)
        frame.pack(fill=tk.BOTH, expand=True)

        data = self.parent.data
        row = 0

        # Bandpass Retention section
        ttk.Label(frame, text="Bandpass Retention:", font=('TkDefaultFont', 10, 'bold')).grid(
            row=row, column=0, columnspan=4, sticky=tk.W, pady=(5, 2))
        row += 1

        ttk.Label(frame, text="Inner:").grid(row=row, column=0, sticky=tk.W, pady=2)
        self.bp_var = tk.StringVar(value='%.2f' % data.bandpassCut)
        self.bp_entry = ttk.Entry(frame, textvariable=self.bp_var, width=10, state='readonly')
        self.bp_entry.grid(row=row, column=1, pady=2, padx=2)
        ttk.Button(frame, text="-", width=3, command=self.onBPDecrease).grid(row=row, column=2, pady=2, padx=2)
        ttk.Button(frame, text="+", width=3, command=self.onBPIncrease).grid(row=row, column=3, pady=2, padx=2)
        row += 1

        # Drift Curve section
        ttk.Label(frame, text="Drift Curve:", font=('TkDefaultFont', 10, 'bold')).grid(
            row=row, column=0, columnspan=4, sticky=tk.W, pady=(10, 2))
        row += 1

        ttk.Label(frame, text="Fit order:").grid(row=row, column=0, sticky=tk.W, pady=2)
        self.drift_order_var = tk.StringVar(value='%i' % data.driftOrder)
        self.drift_order_entry = ttk.Entry(frame, textvariable=self.drift_order_var, width=10, state='readonly')
        self.drift_order_entry.grid(row=row, column=1, pady=2, padx=2)
        ttk.Button(frame, text="-", width=3, command=self.onDriftOrderDecrease).grid(row=row, column=2, pady=2, padx=2)
        ttk.Button(frame, text="+", width=3, command=self.onDriftOrderIncrease).grid(row=row, column=3, pady=2, padx=2)
        row += 1

        ttk.Label(frame, text="Threshold:").grid(row=row, column=0, sticky=tk.W, pady=2)
        self.drift_cut_var = tk.StringVar(value='%i' % data.driftCut)
        self.drift_cut_entry = ttk.Entry(frame, textvariable=self.drift_cut_var, width=10, state='readonly')
        self.drift_cut_entry.grid(row=row, column=1, pady=2, padx=2)
        ttk.Button(frame, text="-", width=3, command=self.onDriftCutDecrease).grid(row=row, column=2, pady=2, padx=2)
        ttk.Button(frame, text="+", width=3, command=self.onDriftCutIncrease).grid(row=row, column=3, pady=2, padx=2)
        row += 1

        # Spectral Kurtosis section
        ttk.Label(frame, text="Spectral Kurtosis:", font=('TkDefaultFont', 10, 'bold')).grid(
            row=row, column=0, columnspan=4, sticky=tk.W, pady=(10, 2))
        row += 1

        ttk.Label(frame, text="Sections:").grid(row=row, column=0, sticky=tk.W, pady=2)
        self.kurt_sec_var = tk.StringVar(value='%i' % data.kurtosisSec)
        self.kurt_sec_entry = ttk.Entry(frame, textvariable=self.kurt_sec_var, width=10, state='readonly')
        self.kurt_sec_entry.grid(row=row, column=1, pady=2, padx=2)
        ttk.Button(frame, text="-", width=3, command=self.onKurtSecDecrease).grid(row=row, column=2, pady=2, padx=2)
        ttk.Button(frame, text="+", width=3, command=self.onKurtSecIncrease).grid(row=row, column=3, pady=2, padx=2)
        row += 1

        ttk.Label(frame, text="Threshold:").grid(row=row, column=0, sticky=tk.W, pady=2)
        self.kurt_cut_var = tk.StringVar(value='%i' % data.kurtosisCut)
        self.kurt_cut_entry = ttk.Entry(frame, textvariable=self.kurt_cut_var, width=10, state='readonly')
        self.kurt_cut_entry.grid(row=row, column=1, pady=2, padx=2)
        ttk.Button(frame, text="-", width=3, command=self.onKurtCutDecrease).grid(row=row, column=2, pady=2, padx=2)
        ttk.Button(frame, text="+", width=3, command=self.onKurtCutIncrease).grid(row=row, column=3, pady=2, padx=2)
        row += 1

        # Separator
        ttk.Separator(frame, orient=tk.HORIZONTAL).grid(row=row, column=0, columnspan=4, sticky='ew', pady=10)
        row += 1

        # Ok button
        ttk.Button(frame, text="Ok", command=self.onOk).grid(row=row, column=3, pady=5)

    def onBPDecrease(self):
        data = self.parent.data
        if data.bandpassCut > 0.05:
            data.bandpassCut -= 0.05
            self.bp_var.set('%.2f' % data.bandpassCut)

    def onBPIncrease(self):
        data = self.parent.data
        if data.bandpassCut < 1.0:
            data.bandpassCut += 0.05
            self.bp_var.set('%.2f' % data.bandpassCut)

    def onDriftOrderDecrease(self):
        data = self.parent.data
        if data.driftOrder > 1:
            data.driftOrder -= 1
            self.drift_order_var.set('%i' % data.driftOrder)

    def onDriftOrderIncrease(self):
        data = self.parent.data
        if data.driftOrder < 12:
            data.driftOrder += 1
            self.drift_order_var.set('%i' % data.driftOrder)

    def onDriftCutDecrease(self):
        data = self.parent.data
        if data.driftCut > 2:
            data.driftCut -= 1
            self.drift_cut_var.set('%i' % data.driftCut)

    def onDriftCutIncrease(self):
        data = self.parent.data
        import math
        max_cut = math.ceil(data.spec.shape[1] / 300)
        if data.driftCut < max_cut:
            data.driftCut += 1
            self.drift_cut_var.set('%i' % data.driftCut)

    def onKurtSecDecrease(self):
        data = self.parent.data
        if data.kurtosisSec > 1:
            data.kurtosisSec -= 1
            self.kurt_sec_var.set('%i' % data.kurtosisSec)

    def onKurtSecIncrease(self):
        data = self.parent.data
        if data.kurtosisSec < 45:
            data.kurtosisSec += 1
            self.kurt_sec_var.set('%i' % data.kurtosisSec)

    def onKurtCutDecrease(self):
        data = self.parent.data
        if data.kurtosisCut > 2:
            data.kurtosisCut -= 1
            self.kurt_cut_var.set('%i' % data.kurtosisCut)

    def onKurtCutIncrease(self):
        data = self.parent.data
        if data.kurtosisCut < 12:
            data.kurtosisCut += 1
            self.kurt_cut_var.set('%i' % data.kurtosisCut)

    def onOk(self):
        self.parent.cAdjust = None
        self.destroy()


class WaterfallDisplay(tk.Toplevel):
    """Zoomable waterfall display window"""

    def __init__(self, parent):
        super().__init__(parent)
        self.parent = parent

        self.title("Zoomable Waterfall")
        self.geometry('800x600')

        self._build_ui()
        self._draw()
        self._connect()

    def _build_ui(self):
        # Status bar at bottom
        self.statusbar = ttk.Label(self, text="", relief=tk.SUNKEN, anchor=tk.W)
        self.statusbar.pack(side=tk.BOTTOM, fill=tk.X)

        frame = ttk.Frame(self)
        frame.pack(fill=tk.BOTH, expand=True)

        self.figure = Figure(figsize=(8, 6))
        self.canvas = FigureCanvasTkAgg(self.figure, master=frame)

        # Create toolbar first so it reserves space at bottom
        self.toolbar = NavigationToolbar2Tk(self.canvas, frame)
        self.toolbar.update()
        self.toolbar.set_message = lambda s: None  # Disable default coordinate display

        # Then pack canvas to fill remaining space
        self.canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)

    def _draw(self):
        data = self.parent.data

        if data.index // (data.spec.shape[0]//2) == 0:
            self.freq = data.freq1
        else:
            self.freq = data.freq2

        if data.bandpass:
            spec = data.specBandpass
            limits = data.limitsBandpass
        else:
            spec = data.spec
            limits = data.limits

        self.figure.clf()
        ax = self.figure.gca()

        if data.usedB:
            m = ax.imshow(to_dB(spec[data.index,:,:]), interpolation='nearest',
                          extent=(self.freq[0]/1e6, self.freq[-1]/1e6, data.time[0], data.time[-1]),
                          origin='lower', cmap=data.cmap,
                          norm=data.norm(limits[data.index][0], limits[data.index][1]))
            cbar = self.figure.colorbar(m)
            cbar.ax.set_ylabel('PSD [arb. dB]')
        else:
            m = ax.imshow(spec[data.index,:,:], interpolation='nearest',
                          extent=(self.freq[0]/1e6, self.freq[-1]/1e6, data.time[0], data.time[-1]),
                          origin='lower', cmap=data.cmap,
                          norm=data.norm(limits[data.index][0], limits[data.index][1]))
            cbar = self.figure.colorbar(m)
            cbar.ax.set_ylabel('PSD [arb. lin.]')

        ax.axis('auto')
        ax.set_xlabel('Frequency [MHz]')
        ax.set_ylabel('Elapsed Time [s]')

        if data.linear:
            tun = data.index // 4 + 1
            ind = data.index % 4
            mapper = {0: 'XX', 1: 'Re(XY)', 2: 'Im(XY)', 3: 'YY'}
            ax.set_title('Tuning %i, %s' % (tun, mapper[ind]))
        else:
            tun = data.index // 4 + 1
            ind = data.index % 4
            mapper = {0: 'I', 1: 'Q', 2: 'U', 3: 'V'}
            ax.set_title('Tuning %i, %s' % (tun, mapper[ind]))

        self.figure.tight_layout()
        self.canvas.draw()

    def _connect(self):
        """Connect to motion events for coordinate display"""
        self.cidmotion = self.figure.canvas.mpl_connect('motion_notify_event', self.on_motion)

    def on_motion(self, event):
        """Update status bar with coordinates"""
        data = self.parent.data

        if event.inaxes:
            clickX = event.xdata
            clickY = event.ydata

            dataX = numpy.where(numpy.abs(clickX - self.freq/1e6) == (numpy.abs(clickX - self.freq/1e6).min()))[0][0]
            dataY = numpy.where(numpy.abs(clickY - data.time) == (numpy.abs(clickY - data.time).min()))[0][0]

            if data.usedB:
                value = to_dB(data.spec[data.index, dataY, dataX])
                self.statusbar.config(text="f=%.4f MHz, t=%.4f s, p=%.2f dB" % (clickX, clickY, value))
            else:
                value = data.spec[data.index, dataY, dataX]
                self.statusbar.config(text="f=%.4f MHz, t=%.4f s, p=%.2f" % (clickX, clickY, value))
        else:
            self.statusbar.config(text="")


class DriftCurveDisplay(tk.Toplevel):
    """Zoomable drift curve display window"""

    def __init__(self, parent):
        super().__init__(parent)
        self.parent = parent

        self.title("Zoomable Drift Curve")
        self.geometry('800x600')

        self._build_ui()
        self._draw()
        self._connect()

        # Set up observer for LST calculation
        self.site = stations.lwa1.get_observer()

    def _build_ui(self):
        # Status bar at bottom
        self.statusbar = ttk.Label(self, text="", relief=tk.SUNKEN, anchor=tk.W)
        self.statusbar.pack(side=tk.BOTTOM, fill=tk.X)

        frame = ttk.Frame(self)
        frame.pack(fill=tk.BOTH, expand=True)

        self.figure = Figure(figsize=(8, 6))
        self.canvas = FigureCanvasTkAgg(self.figure, master=frame)

        # Create toolbar first so it reserves space at bottom
        self.toolbar = NavigationToolbar2Tk(self.canvas, frame)
        self.toolbar.update()
        self.toolbar.set_message = lambda s: None  # Disable default coordinate display

        # Then pack canvas to fill remaining space
        self.canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)

    def _draw(self):
        data = self.parent.data

        if data.bandpass:
            spec = data.specBandpass
        else:
            spec = data.spec

        self.drift = spec[:,:,spec.shape[2]//8:7*spec.shape[2]//8].mean(axis=2)

        self.figure.clf()
        ax = self.figure.gca()

        if data.usedB:
            z = to_dB(self.drift[data.index,:])
            ax.scatter(data.time, z, c=z, marker='x', cmap=data.cmap)
            ax.set_ylabel('Inner 75% Mean Power [arb. dB]')
        else:
            z = self.drift[data.index,:]
            ax.scatter(data.time, z, c=z, marker='x', cmap=data.cmap)
            ax.set_ylabel('Inner 75% Mean Power [arb. lin.]')

        # Add line segments connecting points
        levels = []
        segments = []
        for i in range(1, z.size):
            levels.append(0.5*(z[i-1]+z[i]))
            segments.append([(data.time[i-1], z[i-1]), (data.time[i], z[i])])
        from matplotlib.collections import LineCollection
        from matplotlib.colors import Normalize
        lc = LineCollection(segments)
        lc.set_array(numpy.array(levels))
        lc.set_norm(Normalize(vmin=z.min(), vmax=z.max()))
        lc.set_cmap(data.cmap)
        ax.add_collection(lc)

        ax.set_xlim((data.time[0], data.time[-1]))
        ax.set_xlabel('Elapsed Time [s]')

        if data.linear:
            tun = data.index // 4 + 1
            ind = data.index % 4
            mapper = {0: 'XX', 1: 'Re(XY)', 2: 'Im(XY)', 3: 'YY'}
            ax.set_title('Tuning %i, %s' % (tun, mapper[ind]))
        else:
            tun = data.index // 4 + 1
            ind = data.index % 4
            mapper = {0: 'I', 1: 'Q', 2: 'U', 3: 'V'}
            ax.set_title('Tuning %i, %s' % (tun, mapper[ind]))

        self.figure.tight_layout()
        self.canvas.draw()

    def _connect(self):
        """Connect to motion events for coordinate display"""
        self.cidmotion = self.figure.canvas.mpl_connect('motion_notify_event', self.on_motion)

    def on_motion(self, event):
        """Update status bar with coordinates"""
        data = self.parent.data

        if event.inaxes:
            clickX = event.xdata

            dataX = numpy.where(numpy.abs(clickX - data.time) == (numpy.abs(clickX - data.time).min()))[0][0]

            ts = datetime.utcfromtimestamp(data.timesNPZRestricted[dataX])
            self.site.date = ts.strftime('%Y/%m/%d %H:%M:%S')
            lst = self.site.sidereal_time()

            if data.usedB:
                value = to_dB(self.drift[data.index, dataX])
                self.statusbar.config(text="t=%s, LST=%s, p=%.2f dB" % (ts, lst, value))
            else:
                value = self.drift[data.index, dataX]
                self.statusbar.config(text="t=%s, LST=%s, p=%.2f" % (ts, lst, value))
        else:
            self.statusbar.config(text="")


class PowerSpectrumDisplay(tk.Toplevel):
    """Zoomable power spectrum display window (FFT of drift curve)"""

    def __init__(self, parent):
        super().__init__(parent)
        self.parent = parent

        self.title("Zoomable Power Spectrum")
        self.geometry('800x600')

        self._build_ui()
        self._draw()
        self._connect()

    def _build_ui(self):
        # Status bar at bottom
        self.statusbar = ttk.Label(self, text="", relief=tk.SUNKEN, anchor=tk.W)
        self.statusbar.pack(side=tk.BOTTOM, fill=tk.X)

        frame = ttk.Frame(self)
        frame.pack(fill=tk.BOTH, expand=True)

        self.figure = Figure(figsize=(8, 6))
        self.canvas = FigureCanvasTkAgg(self.figure, master=frame)

        # Create toolbar first so it reserves space at bottom
        self.toolbar = NavigationToolbar2Tk(self.canvas, frame)
        self.toolbar.update()
        self.toolbar.set_message = lambda s: None  # Disable default coordinate display

        # Then pack canvas to fill remaining space
        self.canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)

    def _draw(self):
        data = self.parent.data

        if data.bandpass:
            spec = data.specBandpass
        else:
            spec = data.spec

        # Compute FFT of drift curve
        drift = spec[:,:,spec.shape[2]//8:7*spec.shape[2]//8].mean(axis=2)
        self.fft = numpy.abs(numpy.fft.fft(drift, axis=1))**2
        self.fft_freq = numpy.fft.fftfreq(self.fft.shape[1],
                                          d=(data.time[1] - data.time[0]))
        self.fft = self.fft[:,:self.fft.shape[1]//2]
        self.fft_freq = self.fft_freq[:self.fft_freq.size//2]
        self.fft_units = 'Hz'
        if self.fft_freq.max() < 1:
            self.fft_freq *= 1000.0
            self.fft_units = 'mHz'
        elif self.fft_freq.max() > 1000:
            self.fft_freq /= 1000.0
            self.fft_units = 'kHz'

        self.figure.clf()
        ax = self.figure.gca()

        z = to_dB(self.fft[data.index,:])
        ax.scatter(self.fft_freq, z, c=z, marker='x', cmap=data.cmap)
        ax.set_ylabel('PS of Inner 75% Mean Power [arb. dB]')

        # Add line segments connecting points
        levels = []
        segments = []
        for i in range(1, z.size):
            levels.append(0.5*(z[i-1]+z[i]))
            segments.append([(self.fft_freq[i-1], z[i-1]), (self.fft_freq[i], z[i])])
        from matplotlib.collections import LineCollection
        from matplotlib.colors import Normalize
        lc = LineCollection(segments)
        lc.set_array(numpy.array(levels))
        lc.set_norm(Normalize(vmin=z.min(), vmax=z.max()))
        lc.set_cmap(data.cmap)
        ax.add_collection(lc)

        ax.set_xlim((self.fft_freq[0], self.fft_freq[-1]))
        ax.set_xlabel('Frequency [%s]' % self.fft_units)

        if data.linear:
            tun = data.index // 4 + 1
            ind = data.index % 4
            mapper = {0: 'XX', 1: 'Re(XY)', 2: 'Im(XY)', 3: 'YY'}
            ax.set_title('Tuning %i, %s' % (tun, mapper[ind]))
        else:
            tun = data.index // 4 + 1
            ind = data.index % 4
            mapper = {0: 'I', 1: 'Q', 2: 'U', 3: 'V'}
            ax.set_title('Tuning %i, %s' % (tun, mapper[ind]))

        self.figure.tight_layout()
        self.canvas.draw()

    def _connect(self):
        """Connect to motion events for coordinate display"""
        self.cidmotion = self.figure.canvas.mpl_connect('motion_notify_event', self.on_motion)

    def on_motion(self, event):
        """Update status bar with coordinates"""
        data = self.parent.data

        if event.inaxes:
            clickX = event.xdata

            dataX = numpy.where(numpy.abs(clickX - self.fft_freq) == (numpy.abs(clickX - self.fft_freq).min()))[0][0]
            fft_freq = self.fft_freq[dataX]

            value = to_dB(self.fft[data.index, dataX])
            self.statusbar.config(text="f=%.3f %s, p=%.2f dB" % (fft_freq, self.fft_units, value))
        else:
            self.statusbar.config(text="")


class HelpWindow(tk.Toplevel):
    """
    Help window displaying HTML documentation.

    Provides basic HTML rendering using tk.Text with tags, supporting:
    - Headers (h4, h6)
    - Bold, italic, underline text
    - Unordered and ordered lists
    - Internal anchor links (scrolls to section)
    - External links (opens in browser)
    """

    HELP_HTML = """<html>
<body>
<a name="top"><h4>Table of Contents</h4></a>
<ul>
<li><a href="#intro">Introduction</a></li>
<li><a href="#layout">Window Layout</a></li>
<li><a href="#usage">Usage</a></li>
<li><a href="#mouse">Mouse Interaction</a></li>
<li><a href="#keyboard">Keyboard Interaction</a></li>
</ul>

<p>
<a name="intro">
<h6>Introduction</h6>
plotHDF is a graphical interface for working with HDF5 files created by the hdfWaterfall or drspec2hdf
scripts.  This script allows the user to look at dynamic waterfalls, examine the power as a function
of time, mask RFI, and fit power laws to the data.
<br /><a href="#top">Top</a>
</a>
</p>

<p>
<a name="layout">
<h6>Window Layout</h6>
The plotHDF window is broken into two vertical sections.  The top section shows, from left to right, the
dynamic waterfall, the X and Y polarization data saturation fraction, and the frequency-integrated power
as a function of time.  Any of these three panels can be clicked with the mouse to select a particular
time section to be detailed in the lower section.<br /><br />
The lower section also allows the users to interact with the spectrum, e.g., shifted, zoomed, etc., using the
standard matplotlib toolbar.
<br /><a href="#top">Top</a>
</a>
</p>

<p>
<a name="usage">
<h6>Usage</h6>
After the HDF5 file has been loaded there are a variety of menu, <a href="#mouse">mouse</a>, and
<a href="#keyboard">keyboard commands</a> that can be used to interact with the data.  The menus are:
<ul>
    <li>Color - Adjust the color stretch, map, and transfer function used in the dynamic waterfall</li>
    <li>Data - Change which tuning and data product is currently displayed</li>
    <li>Mask - Use (pseudo) spectral kurtosis to generate a radio frequency interference (RFI) mask and
    apply it to the data.  The masking parameters can also be tweaked in this menu.</li>
    <li>Bandpass - Apply a data or instrumental spectral bandpass to the data.</li>
    <li>Details - Get metadata about the current HDF5 file and how the data were reduced and display
    zoom-able dynamic waterfall, drift curve, and drift curve power spectrum plots.</li>
    <li>Help - Show this help message.</li>
</ul>
<br /><br />
<b>Note:</b> The RFI masks generated either through the "Mask" menu or the mouse/keyboard can be saved
to the HDF5 file at any time using the "File" menu.
<br /><a href="#top">Top</a>
</a>
</p>

<p>
<a name="mouse">
<h6>Mouse Interaction</h6>
The mouse can be used in any of the plotting panels.  In the upper section:
<ul>
    <li>Left Click - Select a particular time to be displayed in the lower spectrum window</li>
    <li>Middle Click - Unmask the timestamp currently under the cursor</li>
    <li>Right Clock - Mask the timestamp currently under the cursor</li>
</ul>
<br /><br />
In the lower section:
<ul>
    <li>Left Click - Select the lower panel for interaction/interact with the matplotlib toolbar</li>
    <li>Middle Click - Unmask the frequency bin currently under the cursor</li>
    <li>Right Click - Mask the frequency bin currently under the cursor</li>
</ul>
<br /><a href="#top">Top</a>
</a>
</p>

<p>
<a name="keyboard">
<h6>Keyboard Interaction</h6>
The keyboard can also be used to interact with any of the plotting panels.  The keyboard commands
are:
<ul>
    <li>Up or Right Arrow - Move the spectrum display forward by one time step <i>(upper section only)</i></li>
    <li>Down or Left Arrow - Move the spectrum display backward by one time step <i>(upper section only)</i></li>
    <li>m - Mask the current timestamp (upper section) or frequency bin (lower section)</li>
    <li>u - Unmask the current timestamp (upper section) or frequency bin (lower section)</li>
    <li>p - Print the power for the current frequency bin <i>(lower section only)</i></li>
    <li>s - Print statistics about the current frequency bin <i>(lower section only)</i></li>
    <li>f - Pick a frequency boundary for a power law fit <i>(lower section only)</i>  <b>Note:</b> This needs to be be done twice to perform the fit</li>
    <li>c - Clear the current power law fit <i>(lower section only)</i></li>
</ul>
<br /><br />
<b>Note:</b> In order to interact with the lower section you will need to click in the axes with the left mouse button.
<br /><a href="#top">Top</a>
</a>
</p>

</body>
</html>"""

    def __init__(self, parent):
        super().__init__(parent)
        self.parent = parent

        self.title("plotHDF Handbook")
        self.geometry('600x500')

        # Store anchor positions for internal navigation
        self.anchors = {}

        self._build_ui()

    def _build_ui(self):
        frame = ttk.Frame(self, padding=5)
        frame.pack(fill=tk.BOTH, expand=True)

        # Scrollable text widget
        text_frame = ttk.Frame(frame)
        text_frame.pack(fill=tk.BOTH, expand=True)

        scrollbar = ttk.Scrollbar(text_frame, orient=tk.VERTICAL)
        scrollbar.pack(side=tk.RIGHT, fill=tk.Y)

        self.text = tk.Text(text_frame, wrap=tk.WORD, padx=10, pady=10,
                           yscrollcommand=scrollbar.set)
        self.text.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        scrollbar.config(command=self.text.yview)

        # Configure text tags for HTML rendering
        self._configure_tags()

        # Render help content
        self._render_html(self.HELP_HTML)
        self.text.config(state='disabled')

    def _configure_tags(self):
        """Configure text tags for HTML-like formatting."""
        import tkinter.font as tkfont

        # Get default font and create variants
        default_font = tkfont.nametofont('TkDefaultFont')
        base_size = default_font.cget('size')
        base_family = default_font.cget('family')

        # Header tags
        self.text.tag_configure('h4', font=(base_family, base_size + 4, 'bold'),
                                spacing1=10, spacing3=5)
        self.text.tag_configure('h6', font=(base_family, base_size + 2, 'bold'),
                                spacing1=8, spacing3=3)

        # Text formatting tags
        self.text.tag_configure('bold', font=(base_family, base_size, 'bold'))
        self.text.tag_configure('italic', font=(base_family, base_size, 'italic'))
        self.text.tag_configure('underline', underline=True)

        # Superscript and subscript (adjust offset)
        self.text.tag_configure('sup', offset=4, font=(base_family, base_size - 2))
        self.text.tag_configure('sub', offset=-4, font=(base_family, base_size - 2))

        # List item tags (with indent)
        self.text.tag_configure('listitem', lmargin1=20, lmargin2=35)
        self.text.tag_configure('olitem', lmargin1=20, lmargin2=35)

        # Link tags
        self.text.tag_configure('link', foreground='blue', underline=True)
        self.text.tag_bind('link', '<Enter>', lambda e: self.text.config(cursor='hand2'))
        self.text.tag_bind('link', '<Leave>', lambda e: self.text.config(cursor=''))

    def _render_html(self, html_content):
        """Parse and render HTML content to the text widget."""
        import re

        # Process the HTML content
        # Normalize multiple whitespace to single space, collapse newlines
        html_content = re.sub(r'\s+', ' ', html_content)

        # Find all tags and text segments
        pattern = r'(<[^>]+>|[^<]+)'
        segments = re.findall(pattern, html_content)

        # Stack to track active formatting
        format_stack = []

        # Link tracking
        current_link_href = None
        link_counter = 0

        # Ordered list tracking
        ol_counter = 0
        in_ol = False

        # Track if we just inserted a newline (to avoid redundant spaces after newlines)
        last_was_newline = True

        for segment in segments:
            if not segment or not segment.strip():
                continue

            if segment.startswith('<'):
                # Strip whitespace from tags only
                segment = segment.strip()
                tag_match = re.match(r'<(/?)(\w+)([^>]*)/?>', segment)
                if not tag_match:
                    continue

                is_closing = tag_match.group(1) == '/'
                tag_name = tag_match.group(2).lower()
                tag_attrs = tag_match.group(3)

                # Handle different tags
                if tag_name in ('html', 'body'):
                    continue

                elif tag_name == 'h4':
                    if not is_closing:
                        format_stack.append('h4')
                    else:
                        if 'h4' in format_stack:
                            format_stack.remove('h4')
                        self.text.insert(tk.END, '\n')
                        last_was_newline = True

                elif tag_name == 'h6':
                    if not is_closing:
                        format_stack.append('h6')
                    else:
                        if 'h6' in format_stack:
                            format_stack.remove('h6')
                        self.text.insert(tk.END, '\n')
                        last_was_newline = True

                elif tag_name == 'b':
                    if not is_closing:
                        format_stack.append('bold')
                    elif 'bold' in format_stack:
                        format_stack.remove('bold')

                elif tag_name == 'i':
                    if not is_closing:
                        format_stack.append('italic')
                    elif 'italic' in format_stack:
                        format_stack.remove('italic')

                elif tag_name == 'u':
                    if not is_closing:
                        format_stack.append('underline')
                    elif 'underline' in format_stack:
                        format_stack.remove('underline')

                elif tag_name == 'sup':
                    if not is_closing:
                        format_stack.append('sup')
                    elif 'sup' in format_stack:
                        format_stack.remove('sup')

                elif tag_name == 'sub':
                    if not is_closing:
                        format_stack.append('sub')
                    elif 'sub' in format_stack:
                        format_stack.remove('sub')

                elif tag_name == 'ul':
                    if not is_closing:
                        self.text.insert(tk.END, '\n')
                        last_was_newline = True
                    else:
                        self.text.insert(tk.END, '\n')
                        last_was_newline = True

                elif tag_name == 'ol':
                    if not is_closing:
                        in_ol = True
                        ol_counter = 0
                        self.text.insert(tk.END, '\n')
                        last_was_newline = True
                    else:
                        in_ol = False
                        self.text.insert(tk.END, '\n')
                        last_was_newline = True

                elif tag_name == 'li':
                    if not is_closing:
                        if in_ol:
                            ol_counter += 1
                            self.text.insert(tk.END, f'\n  {ol_counter}. ', 'olitem')
                        else:
                            self.text.insert(tk.END, '\n  \u2022 ', 'listitem')
                        last_was_newline = False  # We have content after the bullet

                elif tag_name == 'br':
                    self.text.insert(tk.END, '\n')
                    last_was_newline = True

                elif tag_name == 'p':
                    if is_closing:
                        self.text.insert(tk.END, '\n\n')
                        last_was_newline = True

                elif tag_name == 'a':
                    if not is_closing:
                        # Extract href and name attributes
                        href_match = re.search(r'href=["\']([^"\']+)["\']', tag_attrs)
                        name_match = re.search(r'name=["\']([^"\']+)["\']', tag_attrs)

                        if name_match:
                            # This is an anchor definition - store current position
                            anchor_name = name_match.group(1)
                            self.anchors[anchor_name] = self.text.index(tk.END)

                        if href_match:
                            # This is a link - add 'link' to format stack and track href
                            current_link_href = href_match.group(1)
                            link_counter += 1
                            # Create unique tag for this specific link
                            link_tag = f'link_{link_counter}'
                            format_stack.append('link')
                            format_stack.append(link_tag)

                            # Pre-configure the click binding for this link tag
                            if current_link_href.startswith('#'):
                                anchor_name = current_link_href[1:]
                                self.text.tag_bind(link_tag, '<Button-1>',
                                    lambda e, a=anchor_name: self._scroll_to_anchor(a))
                            else:
                                url = current_link_href
                                self.text.tag_bind(link_tag, '<Button-1>',
                                    lambda e, u=url: webbrowser.open(u))
                    else:
                        # Remove link tags from format stack
                        if current_link_href:
                            # Remove the unique link tag
                            for tag in list(format_stack):
                                if tag.startswith('link_'):
                                    format_stack.remove(tag)
                                    break
                            # Remove 'link' tag
                            if 'link' in format_stack:
                                format_stack.remove('link')
                            current_link_href = None

            else:
                # This is text content - decode HTML entities
                text = segment
                text = text.replace('&quot;', '"')
                text = text.replace('&amp;', '&')
                text = text.replace('&lt;', '<')
                text = text.replace('&gt;', '>')
                text = text.replace('&nbsp;', ' ')

                if text.strip():
                    # Handle leading space - preserve it unless we just had a newline
                    if text.startswith(' ') and not last_was_newline:
                        leading_space = ' '
                        text = text.lstrip()
                    else:
                        leading_space = ''
                        text = text.lstrip()

                    # Preserve trailing space
                    trailing_space = ' ' if segment.endswith(' ') else ''
                    text = text.rstrip()

                    # Apply current formatting tags
                    tags = tuple(format_stack) if format_stack else ()
                    if leading_space:
                        self.text.insert(tk.END, leading_space)
                    self.text.insert(tk.END, text, tags)
                    if trailing_space:
                        self.text.insert(tk.END, trailing_space)

                    last_was_newline = False

    def _scroll_to_anchor(self, anchor_name):
        """Scroll to an internal anchor."""
        if anchor_name in self.anchors:
            self.text.see(self.anchors[anchor_name])


def main():
    """Main entry point"""
    parser = argparse.ArgumentParser(
        description='View and edit DRX HDF5 waterfall files',
        epilog='NOTE: The bandpass provided by the -I/--instrumental flag is based on '
               'the current "best knowledge" of LWA1 but may not remove all bandpass-related '
               'features in the data.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('filename', type=str, nargs='?',
                        help='HDF5 file to open')
    parser.add_argument('-s', '--skip', type=aph.positive_or_zero_float, default=0.0,
                        help='skip period in seconds before displaying')
    parser.add_argument('-d', '--duration', type=aph.positive_float, default=-1,
                        help='number of seconds to display')
    parser.add_argument('-o', '--observation', type=aph.positive_int, default=1,
                        help='plot the specified observation')
    parser.add_argument('-I', '--instrumental', action='store_true',
                        help='use an instrument-based bandpass composed of the antenna '
                             'impedance mis-match, the ARX response, and the DRX filter '
                             'coefficients')
    fgroup = parser.add_mutually_exclusive_group(required=False)
    fgroup.add_argument('-n', '--split', action='store_true', default=True,
                        help='take ARX to be in the split bandwidth setting for the '
                             'instrument-based bandpass')
    fgroup.add_argument('-f', '--full', action='store_true',
                        help='take ARX to be in the full bandwidth setting for the '
                             'instrument-based bandpass')
    fgroup.add_argument('-r', '--reduced', action='store_true',
                        help='take ARX to be in the reduced bandwidth setting for the '
                             'instrument-based bandpass')
    args = parser.parse_args()

    # Determine bandpass type
    bandpassType = 'data'
    if args.instrumental:
        bandpassType = 'instrumental'

    # Determine ARX filter
    arxFilter = 'split'
    if args.full:
        arxFilter = 'full'
    elif args.reduced:
        arxFilter = 'reduced'

    app = MainWindow()
    app.offset = args.skip
    app.duration = args.duration

    app.render()

    # Create data handler
    app.data = Waterfall_GUI(app, bandpassType=bandpassType, arxFilter=arxFilter)

    # Disable menus until file is loaded
    app._enableMenus(False)
    app.setSaveButton()

    # Load file if provided
    if args.filename:
        app.dirname = os.path.dirname(args.filename)
        app.filename = args.filename

        app.config(cursor='watch')
        app.showLoading("Loading " + os.path.basename(args.filename))

        try:
            app.data.loadData(args.filename, obsID=args.observation)
            app.data.render()
            app.data.draw()

            app._enableMenus(True)
            app.setDataMenuOptions()
        except Exception as e:
            messagebox.showerror("Error", f"Failed to load file: {e}")

        app.hideLoading()
        app.config(cursor='')

    app.mainloop()


if __name__ == '__main__':
    main()
