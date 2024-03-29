#!/usr/bin/env python3

"""
Given the NPZ file generated by burstMovie.py (as opposed to the movie), perform
multilateration (http://en.wikipedia.org/wiki/Multilateration) to locate the 
source of the RFI.
"""

import sys
import numpy

from lsl.common import stations
from lsl.common.dp import fS
from lsl.reader import tbw
from lsl.reader import errors
from lsl.correlator import fx as fxc
from lsl.astro import unix_to_utcjd, DJD_OFFSET

from matplotlib import pyplot as plt
from scipy.optimize import leastsq, fmin


def crossCorrelate(sig1, sig2):
    """
    Cross-correlate two signals to get the lag between the two
    in samples.  Returns a two-element tuple of the lag values 
    in samples and the strength of the correlation.
    """
    
    cc = numpy.fft.fft(sig1)*numpy.fft.fft(sig2).conj()
    cc = numpy.fft.fftshift(numpy.fft.ifft(cc).real)
    lag = -numpy.arange(-len(cc)/2,len(cc)/2)
    
    return lag, cc


def main(args):
    dataDict = numpy.load(args[0])
    data = dataDict['data']
    ssmif = dataDict['ssmif']
    
    if ssmif != '':
        station = stations.parse_ssmif(ssmif)
    else:
        station = stations.lwa1
    antennas = station.antennas

    # Find the time of arrival for all of the various pulses relative
    # to ccPoint
    p = []
    for ccPoint in range(260):
        delays = numpy.zeros(data.shape[0]) - 100
        g = []
        for i in range(delays.size//2):
            lagX, ccX = crossCorrelate(data[2*ccPoint+0,:], data[2*i+0,:])
            lagY, ccY = crossCorrelate(data[2*ccPoint+1,:], data[2*i+1,:])
        
            lX = lagX[numpy.where( ccX == ccX.max() )[0][0]]
            lY = lagY[numpy.where( ccY == ccY.max() )[0][0]]
        
            mX = numpy.log10(ccX.max())
            mY = numpy.log10(ccY.max())
            #print(lX, mX)
        
            if i == 9:
                continue
            #if antennas[2*i+0].stand.id > 256:
                #continue
            
            if mX > 5.7:
                delays[2*i+0] = lagX[numpy.where( ccX == ccX.max() )[0][0]] / fS
            if mY > 5.7:
                delays[2*i+1] = lagY[numpy.where( ccY == ccY.max() )[0][0]] / fS
            
            #print(2*i+0, antennas[2*i+0].stand.id, delays[2*i+0]*1e9, mX)
            #print(2*i+1, antennas[2*i+1].stand.id, delays[2*i+1]*1e9, mY)

            g.append(mY)
        g = numpy.array(g)
        p.append(g.mean())
    p = numpy.array(p)
    print(p.argmax())
        
    ccPoint = p.argmax()
    delays = numpy.zeros(data.shape[0]) - 100
    for i in range(delays.size//2):
        lagX, ccX = crossCorrelate(data[2*ccPoint+0,:], data[2*i+0,:])
        lagY, ccY = crossCorrelate(data[2*ccPoint+1,:], data[2*i+1,:])
    
        lX = lagX[numpy.where( ccX == ccX.max() )[0][0]]
        lY = lagY[numpy.where( ccY == ccY.max() )[0][0]]
    
        mX = numpy.log10(ccX.max())
        mY = numpy.log10(ccY.max())
        print(lX, mX)
    
        if i == 9:
            continue
        #if antennas[2*i+0].stand.id > 256:
            #continue
        
        if mX > 5.7:
            delays[2*i+0] = lagX[numpy.where( ccX == ccX.max() )[0][0]] / fS
        if mY > 5.7:
            delays[2*i+1] = lagY[numpy.where( ccY == ccY.max() )[0][0]] / fS
        
        print(2*i+0, antennas[2*i+0].stand.id, delays[2*i+0]*1e9, mX)
        print(2*i+1, antennas[2*i+1].stand.id, delays[2*i+1]*1e9, mY)

    # Create arrays for the stand ID numbers (ids), combined status codes 
    # (status), polarizations (pols), and stand positions (standPos)
    ids = numpy.array([a.stand.id for a in antennas])
    status = numpy.array([a.combined_status for a in antennas])
    pols = numpy.array([a.pol for a in antennas])
    standPos = numpy.array([[ant.stand.x, ant.stand.y, ant.stand.z] for ant in antennas])

    # Create the master output figure
    fig = plt.figure(figsize=(10,10))
    ax1 = fig.add_subplot(2, 2, 1)
    ax2 = fig.add_subplot(2, 2, 2)
    ax3 = fig.add_subplot(2, 3, 4)
    ax4 = fig.add_subplot(2, 3, 5)
    ax5 = fig.add_subplot(2, 3, 6)
    good = numpy.where( (delays > -1) & (pols == 0) & numpy.isfinite(delays) )[0]
    CB = ax1.scatter(standPos[good,0], standPos[good,1], c=delays[good]*1e9)
    c = fig.colorbar(CB, ax=ax1)
    c.ax.set_ylabel('Time of Arrival [ns]')
    good = numpy.where( (delays > -1) & (pols == 1) & numpy.isfinite(delays) )[0]
    CB = ax2.scatter(standPos[good,0], standPos[good,1], c=delays[good]*1e9)
    c = fig.colorbar(CB, ax=ax2)
    c.ax.set_ylabel('Time of Arrival [ns]')
    
    for ax in [ax1, ax2]:
        ax.plot(0, 0, marker='+',color='k')
        ax.plot([-59.827, 59.771, 60.148, -59.700, -59.827], 
                [59.752, 59.864, -59.618, -59.948, 59.752], linestyle='--', color='k')
        ax.plot([55.863, 58.144, 58.062, 55.791, 55.863], 
                [45.946, 45.999, 51.849, 51.838, 45.946], linestyle='-', color='k')
                
    ax1.set_title('X pol.')
    ax2.set_title('Y pol.')

    # To help with convergence with the limit range in z, zero the z values 
    # for all of the stands.
    standPos[:,2] *= 0

    print(standPos.shape, delays.shape, pols.shape)
    
    analysisRange = 200

    for p in [0,1]:
        valid = numpy.where( (delays > -100) & (pols == p) & (ids <= 256) )[0]
        toUse = numpy.where( (delays > -100) & (pols == p) )[0]
        tau = delays - delays[ccPoint]

        relPos = standPos - standPos[ccPoint,:]

        v = 100.0 / (delays[valid].max() - delays[valid].min())
        #v = 2.9e8
        print(v, v/3e8)
        A = 2.0/v*(relPos[toUse[2:],0] / tau[toUse[2:]] - relPos[toUse[1],0] / tau[toUse[1]])
        B = 2.0/v*(relPos[toUse[2:],1] / tau[toUse[2:]] - relPos[toUse[1],1] / tau[toUse[1]])
        C = 2.0/v*(relPos[toUse[2:],2] / tau[toUse[2:]] - relPos[toUse[1],2] / tau[toUse[1]])
        D =  v*tau[toUse[2:]] - \
            v*tau[toUse[1]] - \
            (relPos[toUse[2:],0]**2 + relPos[toUse[2:],1]**2 + relPos[toUse[2:],2]**2) / v/tau[toUse[2:]] + \
            (relPos[toUse[1], 0]**2 + relPos[toUse[1], 1]**2 + relPos[toUse[1], 2]**2) / v/tau[toUse[1]]

        def toa(pos, values):
            good = numpy.where( numpy.isfinite(values[0,:]) )[0]
            return pos[0]*values[0,good] + pos[1]*values[1,good] + pos[2]*values[2,good] + values[3,good]
            
        def errorFunc(pos, values, D):
            return toa(pos, values)
            
        def errorFunc2(pos, values, D):
            return (toa(pos, values)**2).sum()


        print("Localizing pol %i using %i antennas" % (p, len(toUse)))
        emitter = [500, 500, 0]
        out = fmin(errorFunc2, emitter, args=(numpy.array([A, B, C, D]), numpy.zeros_like(A)), maxiter=100000, maxfun=100000)
        print(out)
        out, cov, ti, tm, te = leastsq(errorFunc, out, args=(numpy.array([A, B, C, D]), numpy.zeros_like(A)), maxfev=100000, full_output=True)
        
        out += standPos[ccPoint,:]
        #print(out[0]*A + out[1]*B + out[2]*C + D)
        try:
            err = numpy.sqrt(cov)
        except:
            err = numpy.zeros((3,3))
        
        print("X: %+.2f +/- %.2f\nY: %+.2f +/- %.2f\nZ: %+.2f +/- %.2f" % (out[0], err[0,0], out[1], err[1,1], out[2], err[2,2]))
        
        x = numpy.linspace(-analysisRange,analysisRange,200)
        y = numpy.linspace(-analysisRange,analysisRange,200)
        Y, X = numpy.meshgrid(x,y)
        out = numpy.zeros_like(X)
        gV = numpy.where( numpy.isfinite(A) )[0]
        for i in range(x.size):
            for j in range(y.size):
                xV = X[i,j] - standPos[ccPoint,0]
                yV = Y[i,j] - standPos[ccPoint,1]
                
                out[i,j] = numpy.log10(((xV*A[gV] + yV*B[gV] + 0*C[gV] + D[gV])**2).sum())*10
                
        try:
            master += 10**(out/10)
        except NameError:
            master = 10**(out/10)
        
        if p == 0:
            ax = ax3
        else:
            ax = ax4
        ax.plot(0,0,marker='+',color='k')
        ax.plot([-59.827, 59.771, 60.148, -59.700, -59.827], 
                [59.752, 59.864, -59.618, -59.948, 59.752], linestyle='--', color='k')
        ax.plot([55.863, 58.144, 58.062, 55.791, 55.863], 
                [45.946, 45.999, 51.849, 51.838, 45.946], linestyle='-', color='k')
        ax.imshow(out, extent=(x.min(), x.max(), y.min(), y.max()), origin='lower')

    ax5.plot(0,0,marker='+',color='k')
    ax5.plot([-59.827, 59.771, 60.148, -59.700, -59.827], 
            [59.752, 59.864, -59.618, -59.948, 59.752], linestyle='--', color='k')
    ax5.plot([55.863, 58.144, 58.062, 55.791, 55.863], 
            [45.946, 45.999, 51.849, 51.838, 45.946], linestyle='-', color='k')
    out = numpy.log10(master)*10
    ax5.imshow(out, extent=(x.min(), x.max(), y.min(), y.max()), origin='lower')

    ax3.set_title('X pol.')
    ax4.set_title('Y pol.')
    ax5.set_title('Combined.')
    
    for ax in [ax1, ax2, ax3, ax4, ax5]:
        ax.set_xlabel('$\Delta$ x [m]')
        ax.set_ylabel('$\Delta$ y [m]')

    plt.show()


if __name__ == "__main__":
    main(sys.argv[1:])
