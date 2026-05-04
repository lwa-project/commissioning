#!/usr/bin/env python3

import os
import sys
import numpy as np
import argparse

from lsl.common.stations import parse_ssmif

from matplotlib import pyplot as plt


def main(args):
    NINPUT_CHUNK = 32
    if args.dig_board:
        NINPUT_CHUNK = 16
        
    delays = []
    for filename in args.filename:
        delays.append([])
        station = parse_ssmif(filename)
        for ant in station.antennas:
            delays[-1].append(ant.cable.delay(74e6))
    delays = np.array(delays)
    nfile, ninput = delays.shape
    
    fig = plt.figure()
    n = len(args.filename)
    for i in range(nfile):
        for j in range(i+1, nfile):
            diff = (delays[j,:] - delays[i,:])*1e9
            
            ax = fig.add_subplot(nfile-1, nfile-1, i*(nfile-1)+j)
            for k in range(0, ninput, NINPUT_CHUNK):
                ax.plot(1+k+np.arange(NINPUT_CHUNK), diff[k:k+NINPUT_CHUNK], linestyle='', marker='x')
            ax.set_xlabel('Input #')
            ax.set_ylabel('Delay Diff [ns]')
    plt.show()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Plot the changes in cable delays between pairs of SSMIFs'
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('filename', type=str, nargs='+',
                        help='SSMIFs to difference')
    parser.add_argument('-d', '--dig-board',  action='store_true',
                        help='color code by digitizer board instead of ZCU102')
    args = parser.parse_args()
    if len(args.filename) < 2:
        raise RuntimeError("Need at least two SSMIFs to difference")
    main(args)
