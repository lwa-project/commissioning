#!/usr/bin/env python3

"""
Read in a SSMIF file and estimate the DRX beam for a given frequency and 
topocentric pointing center.  The estimate is based off a simple delay-and-sum 
beam former so it won't be an exact match.
"""

# Python3 compatiability
from __future__ import print_function, division
import sys
if sys.version_info > (3,):
    xrange = range
    
import os
import sys
import aipy
import time
import numpy
import argparse

from astropy.constants import c as speedOfLight
vLight = speedOfLight.to('m/s').value

from lsl.common.dp import fS
from lsl.common.stations import parse_ssmif
from lsl.misc import beamformer
from lsl.common.paths import DATA as dataPath
from lsl.misc import parser as aph

from matplotlib import pyplot as plt


def main(args):
    station = parse_ssmif(args.filename)
    antennas = station.antennas

    digs    = numpy.array([ant.digitizer  for ant in antennas])
    ants    = numpy.array([ant.id         for ant in antennas])
    stands  = numpy.array([ant.stand.id   for ant in antennas])
    pols    = numpy.array([ant.pol        for ant in antennas])
    antStat = numpy.array([ant.status     for ant in antennas])
    feeStat = numpy.array([ant.fee.status for ant in antennas])

    badStands = numpy.where( antStat != 3 )[0]
    badFees   = numpy.where( feeStat != 3 )[0]
    bad = numpy.where( (stands > 256) | (antStat != 3) | (feeStat != 3) )[0]
    print("Number of bad stands:   %3i" % len(badStands))
    print("Number of bad FEEs:     %3i" % len(badFees))
    print("---------------------------")
    print("Total number bad inputs: %3i" % len(bad))
    print(" ")
    
    # Calculate the beamformer delays in this direction, making sure to
    # do it at 74 MHz and then quantize them.
    delays =  beamformer.calc_delay(antennas, freq=74e6, azimuth=args.azimuth, elevation=args.elevation)
    delays_int = (delays*fS).astype(numpy.int32)
    delays_fra = (16*(delays*fS - delays_int)).astype(numpy.int16)
    delays_quantized = delays_int / fS + delays_fra / 16.0 / fS

    # Build up an array of antenna positions (xyz) and cable delays (cbl)/
    # attenuation (atn) for the frequency of interest.  These will be used 
    # later to compute the physical delays through the system
    xyz = numpy.array([(a.stand.x, a.stand.y, a.stand.z) for a in antennas]).T
    cbl = numpy.array([a.cable.delay((args.frequency*1e6)) for a in antennas])
    atn = numpy.array([a.cable.gain((args.frequency*1e6)) for a in antennas])
    atn = numpy.sqrt(atn / atn.max())
    
    # Apply a uniform weight across the dipoles
    wgt = numpy.ones(len(antennas))

    # Build up lists containing indecies of good antennas in both polarizations
    # and then zero out the weights of bad antennas
    X = [i for i,a in enumerate(antennas) if a.combined_status == 33 and a.pol == 0]
    Y = [i for i,a in enumerate(antennas) if a.combined_status == 33 and a.pol == 1]
    wgt[[i for i,a in enumerate(antennas) if a.combined_status != 33]] = 0.0

    # Quantize the weights since we may need to do that anyways when we throw 
    # this on the station
    wgt_int = (wgt*32767).astype(numpy.int16)
    wgt_quantized = wgt_int / 32767.0

    # Create and azimuth, elevation, and output power arrays
    ires = int(round(1.0/min([1.0, args.resolution])))
    az = numpy.arange(0, 360*ires+1, 1)/float(ires)
    el = numpy.arange(0, 90*ires+1, 1)/float(ires)
    el, az = numpy.meshgrid(el, az)
    pwrX = az*0.0
    pwrY = az*0.0
    
    # Go!
    print("Calculating NS/EW beams for az. %.2f, el. %.2f at %.2f MHz - with %i antennas" % (args.azimuth, args.elevation, args.frequency, wgt.sum()))
    tStart = time.time()
    for i in xrange(az.shape[0]):
        for j in xrange(az.shape[1]):
            ## Get the directon we are intersted in knowing the beam patter in
            a, e = az[i,j]*numpy.pi/180, el[i,j]*numpy.pi/180
            
            ## Convert this direction to a vector and then calculate the physical
            ## delay for a signal from this direction
            pc = numpy.array([numpy.cos(e)*numpy.sin(a), 
                              numpy.cos(e)*numpy.cos(a), 
                              numpy.sin(e)])
            delays_physical = cbl - numpy.dot(pc, xyz) / vLight
            
            ## Calculate the beamformed signal without summing across antennas
            sig = wgt_quantized*atn*numpy.exp(-2j*numpy.pi*(args.frequency*1e6)*(delays_physical-delays_quantized))
            ## Now sum across the right antennas
            pwrX[i,j] = numpy.abs( numpy.sum(sig[X]) )**2
            pwrY[i,j] = numpy.abs( numpy.sum(sig[Y]) )**2
            
    # Calculate the dipole gain pattern to apply as a correction to the beam pattern
    ## Load in the data
    dd = numpy.load(os.path.join(dataPath, 'beam-shape.npz'))
    coeffs = dd['coeffs']
    try:
        dd.close()
    except AttributeError:
        pass
    ## Calculate how many harmonics are stored in the data set and reorder the data
    ## to AIPY's liking
    deg = coeffs.shape[0]-1
    lmax = int((numpy.sqrt(1+8*coeffs.shape[1])-3)/2)
    beamShapeDict = {}
    for i in range(deg+1):
        beamShapeDict[i] = numpy.squeeze(coeffs[-1-i,:])
    ## Build the model
    dipole = aipy.amp.BeamAlm(numpy.array([args.frequency/1e3]), lmax=lmax, mmax=lmax, deg=deg, nside=128, coeffs=beamShapeDict)
    ## Calculate the response for X pol.
    dplX = dipole.response(aipy.coord.azalt2top(numpy.concatenate([[az.ravel()*numpy.pi/180], 
                                                                   [el.ravel()*numpy.pi/180]])))
    ## Rotate the azimuth values by 90 degrees to get Y pol.
    dplY = dipole.response(aipy.coord.azalt2top(numpy.concatenate([[az.ravel()*numpy.pi/180+numpy.pi/2], 
                                                                   [el.ravel()*numpy.pi/180]])))
    ## Re-shape to get back to a 2-D array that matches az (and pwrX/Y)
    dplX.shape = az.shape
    dplY.shape = az.shape
    ## Apply to the beam pattern
    pwrX *= dplX
    pwrY *= dplY
    
    # Save
    for pol,beam in zip(('NS','EW'), (pwrX,pwrY)):
        numpy.savez('%s_%iMHz_%iaz_%iel_%s.npz' % (station.name, args.frequency, args.azimuth, args.elevation, pol), 
                    station=station.name.lower(), beam=beam, freq=(args.frequency*1e6), pol=pol, 
                    az=args.azimuth, el=args.elevation, res=args.resolution)
    print("-> Finished in %.3f seconds" % (time.time() - tStart))
        
    if not args.no_plots:
        norm = max([pwrX.max(), pwrY.max()])
        pwrX /= norm
        pwrY /= norm
        vmin, vmax = 0.0, 1.0
        if args.dB:
            pwrX = numpy.log10(pwrX)*10.0
            pwrY = numpy.log10(pwrY)*10.0
            vmin, vmax = -30.0, 0.0
        fig = plt.figure()
        ax1 = fig.add_subplot(2, 1, 1)
        ax2 = fig.add_subplot(2, 1, 2)
        c = ax1.imshow(pwrX.T, origin='lower', interpolation='nearest', extent=(az.min(), az.max(), el.min(), el.max()), 
                       vmin=vmin, vmax=vmax)
        cb = fig.colorbar(c, ax=ax1)
        cb.set_label('Power'+(' [dB]' if args.dB else ' [lin]'))
        c = ax2.imshow(pwrY.T, origin='lower', interpolation='nearest', extent=(az.min(), az.max(), el.min(), el.max()),
                       vmin=vmin, vmax=vmax)
        cb = fig.colorbar(c, ax=ax2)
        cb.set_label('Power'+(' [dB]' if args.dB else ' [lin]'))
        for ax in (ax1, ax2):
            ax.set_xlabel('Az. [deg]')
            ax.set_ylabel('El. [deg]')
        plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='read in a SSMIF file and estimate the DRX beam for a given frequency and topocentric pointing center', 
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
    parser.add_argument('filename', type=str, 
                        help='SSMIF filename to use')
    parser.add_argument('-f', '--frequency', type=aph.positive_float, default=65.0, 
                        help='frequency in MHz to calculate the beam for')
    parser.add_argument('-a', '--azimuth', type=float, default=90.0, 
                        help='azimuth east of north in degrees for the pointing center')
    parser.add_argument('-e', '--elevation', type=aph.positive_float, default=90.0, 
                        help='elevation above the horizon in degrees for the pointing center')
    parser.add_argument('-r', '--resolution', type=aph.positive_float, default=1.0, 
                        help='beam model resolution in degrees')
    parser.add_argument('-n', '--no-plots', action='store_true', 
                        help='do not show plots of the beam response')
    parser.add_argument('-d', '--dB', action='store_true', 
                        help='show plots with a log stretch')
    args = parser.parse_args()
    main(args)
    
