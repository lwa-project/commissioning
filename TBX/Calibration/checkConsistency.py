#!/usr/bin/env python3

import os
import sys
import numpy as np
import argparse

from lsl.imaging import utils
from lsl.astro import DJD_OFFSET
from lsl.statistics import robust
from lsl.misc import parser as aph

from matplotlib import pyplot as plt


def parse_sc_header(filename):
    """
    Function to parse the text header from a .sc file and return a dictionary
    of its contents.
    """
    
    header = {}
    
    in_settings = False
    with open(filename, 'r') as fh:
        for line in fh:
            line = line.strip()
            if line.find('Settings:') != -1:
                in_settings = True
                continue
            elif line.find('Columns:') != -1:
                in_settings = False
                break
                
            if line.startswith('#') and in_settings:
                try:
                    key, value = line.split(':', 1)
                    key = key.split(None, 1)[-1]
                    value = value.rsplit(None, 1)[0]
                    value = value.strip()
                    if value == 'True':
                        value = True
                    elif value == 'False':
                        value = False
                    else:
                        try:
                            value = int(value, 10)
                        except ValueError:
                            try:
                                value = float(value)
                            except ValueError:
                                pass
                    header[key] = value
                except ValueError as e:
                    pass
                    
    return header


def load_sc_file(filename):
    """
    Wrapper around parse_sc_header() and np.loadtxt() to get everything out
    of a .sc file in a single call.
    """
    
    hdr = parse_sc_header(filename)
    data = np.loadtxt(filename)
    return hdr, data


def main(args):
    # Get the list of .cs filenames to parse
    filenames = args.filename
    
    # Go!
    utcs = []
    lsts = []
    delaysX = {}
    delaysY = {}
    order = []
    n_used = 0
    for filename in filenames:
        ## Figure out what the corresponding FITS IDI file is so that we can pull 
        ## out the date and LST
        fitsname = filename.replace('.sc', '.FITS_1')
        if not os.path.exists(fitsname):
            fitsname = filename.replace('.sc', '.FITS_CAL_1')
        if not os.path.exists(fitsname):
            fitsname = filename.replace('.sc', '.ms_1')
        idi = utils.CorrelatedData(fitsname)
        lo = idi.get_observer()
        lo.date = idi.date_obs.strftime("%Y/%m/%d %H:%M:%S")
        lst = float(lo.sidereal_time()) * 12.0/np.pi
        
        ## Load in the actual file
        hdr, data = load_sc_file(filename)
        if args.only_converged:
            if not (hdr['Converged XX'] and hdr['Converged YY']):
                continue
        utcs.append(lo.date + DJD_OFFSET)
        lsts.append(lst)
        n_used += 1
        
        for i in range(data.shape[0]):
            stand, ax, dx, ay, dy = data[i,:]
            stand = int(stand)
            dx = float(dx)
            dy = float(dy)
            if stand not in order:
                order.append(stand)
                
            try:
                delaysX[stand].append(dx)
                delaysY[stand].append(dy)
            except KeyError:
                delaysX[stand] = [dx,]
                delaysY[stand] = [dy,]
                
    # Convert to arrays
    for stand in delaysX.keys():
        delaysX[stand] = np.array(delaysX[stand])
        valid = np.where( np.isfinite(delaysX[stand]) & (np.abs(delaysX[stand]) <= args.max_delay) )[0]
        if len(valid) < len(delaysX[stand])/2:
            print(stand, 'X')
        if len(valid) > 0:
            try:
                repl = robust.mean(delaysX[stand][valid])
            except (ValueError, ZeroDivisionError):
                repl = np.mean(delaysX[stand][valid])
            delaysX[stand][np.where(~np.isfinite(delaysX[stand]) | (np.abs(delaysX[stand]) > args.max_delay))] = repl
        else:
            delaysX[stand][:] = 0.0
            
        delaysY[stand] = np.array(delaysY[stand])
        valid = np.where( np.isfinite(delaysY[stand]) & (np.abs(delaysY[stand]) <= args.max_delay) )[0]
        if len(valid) < len(delaysX[stand])/2:
            print(stand, 'Y')
        if len(valid) > 0:
            try:
                repl = robust.mean(delaysY[stand][valid])
            except (ValueError, ZeroDivisionError):
                repl = np.mean(delaysY[stand][valid])
            delaysY[stand][np.where(~np.isfinite(delaysY[stand]) | (np.abs(delaysY[stand]) > args.max_delay))] = repl
        else:
            delaysY[stand][:] = 0.0
            
        
    # Calculate the mean delay for each capture
    vsX = np.zeros((len(delaysX.keys()), n_used))
    vsY = np.zeros_like(vsX)
    for i,stand in enumerate(delaysX.keys()):
        dx = delaysX[stand]
        dy = delaysY[stand]
        for j in range(len(dx)):
            vsX[i,j] = dx[j]
            vsY[i,j] = dy[j]
    msX = robust.mean(vsX, axis=0)
    msY = robust.mean(vsY, axis=0)
    
    # Merge the delays together
    with open(args.output, 'w') as fh:
        fh.write("################################\n")
        fh.write("#                              #\n")
        fh.write("# Settings:                    #\n")
        fh.write(f"#  File Count: {len(filenames):4d}            #\n")
        fh.write(f"#  Max Valid Delay: {args.max_delay:6.2f} ns #\n")
        fh.write(f"#  Converged Only: {str(args.only_converged):5s}       #\n")
        fh.write("#                              #\n")
        fh.write("################################\n")
        fh.write("#                              #\n")
        fh.write("# Columns:                     #\n")
        fh.write("# 1) Stand number              #\n")
        fh.write("# 2) Zero                      #\n")
        fh.write("# 3) X pol. delay (ns)         #\n")
        fh.write("# 4) X pol. delay RMS (ns)     #\n")
        fh.write("# 5) Zero                      #\n")
        fh.write("# 6) Y pol. delay (ns)         #\n")
        fh.write("# 7) Y pol. delay RMS (ns)     #\n")
        fh.write("#                              #\n")
        fh.write("################################\n")
        for stand in order:
            try:
                delayX = robust.mean( delaysX[stand] - msX)
                dstdX  = robust.std(  delaysX[stand] - msX)
            except:
                delayX = np.mean( delaysX[stand] - msX)
                dstdX  = np.std(  delaysX[stand] - msX)
            try:
                delayY = robust.mean( delaysY[stand] - msY)
                dstdY  = robust.std(  delaysY[stand] - msY)
            except:
                delayY = np.mean( delaysY[stand] - msY)
                dstdY  = np.std(  delaysY[stand] - msY)
            fh.write("%3i  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f\n" % (stand, 0.0, delayX, dstdX, 0.0, delayY, dstdY))
            
    if args.plot:
        #
        # By LST
        #
        lsts = np.array(lsts)
        orderL = np.argsort(lsts)
        
        figL = plt.figure()
        axLX1 = figL.add_subplot(3, 2, 1)
        axLY1 = figL.add_subplot(3, 2, 2)
        axLX2 = figL.add_subplot(3, 2, 3)
        axLY2 = figL.add_subplot(3, 2, 4)
        axLX3 = figL.add_subplot(3, 2, 5)
        axLY3 = figL.add_subplot(3, 2, 6)
        for stand in delaysX.keys():
            axLX1.plot(lsts[orderL], delaysX[stand][orderL])
            axLY1.plot(lsts[orderL], delaysY[stand][orderL])
        axLX1.set_title('X pol.')
        axLX1.set_xlabel("LST [hours]")
        axLX1.set_ylabel("$\\tau_X$ [ns]")
        axLY1.set_title('Y pol.')
        axLY1.set_xlabel("LST [hours]")
        axLY1.set_ylabel("$\\tau_Y$ [ns]")
        
        axLX2.plot(lsts[orderL], msX[orderL], linestyle='-', marker='x')
        axLY2.plot(lsts[orderL], msY[orderL], linestyle='-', marker='x')
        axLX2.set_xlabel("LST [hours]")
        axLX2.set_ylabel("$|\\tau_x|$ [ns]")
        axLY2.set_xlabel("LST [hours]")
        axLY2.set_ylabel("$|\\tau_y|$ [ns]")
        
        for stand in delaysX.keys():
            dx = delaysX[stand] - msX
            dy = delaysY[stand] - msY
            axLX3.plot(lsts[orderL], dx[orderL])
            axLY3.plot(lsts[orderL], dy[orderL])
        axLX3.set_xlabel("LST [hours]")
        axLX3.set_ylabel("$\\tau_X-|\\tau_x|$ [ns]")
        axLY3.set_xlabel("LST [hours]")
        axLY3.set_ylabel("$\\tau_Y-|\\tau_y|$ [ns]")
        
        figL.tight_layout()
        
        #
        # By JD
        #
        utcs = np.array(utcs)
        orderJ = np.argsort(utcs)
        utcOffset = utcs.min()
        utcs -= utcOffset
        
        figJ = plt.figure()
        axJX1 = figJ.add_subplot(3, 2, 1)
        axJY1 = figJ.add_subplot(3, 2, 2)
        axJX2 = figJ.add_subplot(3, 2, 3)
        axJY2 = figJ.add_subplot(3, 2, 4)
        axJX3 = figJ.add_subplot(3, 2, 5)
        axJY3 = figJ.add_subplot(3, 2, 6)
        for stand in delaysX.keys():
            axJX1.plot(utcs[orderJ], delaysX[stand][orderJ])
            axJY1.plot(utcs[orderJ], delaysY[stand][orderJ])
        axJX1.set_title('X pol.')
        axJX1.set_xlabel("JD [days - %.1f]" % utcOffset)
        axJX1.set_ylabel("$\\tau_X$ [ns]")
        axJY1.set_title('Y pol.')
        axJY1.set_xlabel("JD [days - %.1f]" % utcOffset)
        axJY1.set_ylabel("$\\tau_Y$ [ns]")
        
        axJX2.plot(utcs[orderJ], msX[orderJ], linestyle='-', marker='x')
        axJY2.plot(utcs[orderJ], msY[orderJ], linestyle='-', marker='x')
        axJX2.set_xlabel("JD [days - %.1f]" % utcOffset)
        axJX2.set_ylabel("$|\\tau_x|$ [ns]")
        axJY2.set_xlabel("JD [days - %.1f]" % utcOffset)
        axJY2.set_ylabel("$|\\tau_y|$ [ns]")
        
        tooHighX = []
        tooHighY = []
        for stand in delaysX.keys():
            dx = delaysX[stand] - msX
            dy = delaysY[stand] - msY
            print("%3i:  %6.3f +/- %5.3f  %6.3f +/- %5.3f" % (stand, dx.mean(), dx.std(), dy.mean(), dy.std()))
            if dx.std() > 1:
                tooHighX.append(stand)
            if dy.std() > 1:
                tooHighY.append(stand)
            axJX3.plot(utcs[orderJ], dx[orderJ])
            axJY3.plot(utcs[orderJ], dy[orderJ])
        axJX3.set_xlabel("JD [days - %.1f]" % utcOffset)
        axJX3.set_ylabel("$\\tau_X-|\\tau_x|$ [ns]")
        axJY3.set_xlabel("JD [days - %.1f]" % utcOffset)
        axJY3.set_ylabel("$\\tau_Y-|\\tau_y|$ [ns]")
        
        print(f"RMS >1 ns:  X={len(tooHighX)}, Y={len(tooHighY)}")
        print(f"X: {tooHighX}, Y: {tooHighY}")
        
        figJ.tight_layout()
        
        plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="given a collection of delay files generated by applySelfCalTBX.py, check the internal consistency of the results and updated SSMIF",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
    parser.add_argument('filename', type=str, nargs='+',
                        help='filename to check')
    parser.add_argument('-c', '--only-converged', action='store_true',
                        help='only consider files where the self cal converged in both pols.')
    parser.add_argument('-d', '--max-delay', type=aph.positive_float, default=1000,
                        help='maximum delay in ns to consider valid')
    parser.add_argument('-p', '--plot', action='store_true',
                        help='show diagnostic plots')
    parser.add_argument('-o', '--output', type=str, default='merged.delays',
                        help='write output to the specified filename')
    args = parser.parse_args()
    main(args)
    
