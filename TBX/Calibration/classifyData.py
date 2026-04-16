#!/usr/bin/env python3

import os
import sys
import math
import ephem
import argparse

from lsl.common.stations import parse_ssmif
from lsl.reader.ldp import TBXFile


_srcs = ["CygA,f|J,19:59:28.30,+40:44:02.0,1"]

# We want Cyg data to be within 1 hour of transit
_cyg_transit_window = [-1, 1]

# We want diffuse data to be between the transits of Tau and Vir
_diffuse_lst_window = [6.5, 10.5]


def main(args):
    # Break out the files we need
    ssmif = args.ssmif
    filenames = args.filename
    
    # Setup the LWA station information
    station = parse_ssmif(ssmif)
    antennas = station.antennas
    
    # Get an observer reader for calculations
    obs = station.get_observer()
    
    # Load in the sources and pull out the correct one
    sources = {}
    for line in _srcs:
        bdy = ephem.readdb(line)
        sources[bdy.name] = bdy
    cyg = sources['CygA']
    
    for filename in filenames:
        ## Find the target azimuth/altitude to use
        idf = TBXFile(filename)
        tStart = idf.get_info('start_time').datetime
        freqs = idf.get_info('freq1')
        idf.close()
        
        obs.date = tStart.strftime("%Y/%m/%d %H:%M:%S")
        cyg.compute(obs)
        lst = obs.sidereal_time() * 12/math.pi
        cyg_ha = cyg.ra*12/math.pi - lst
        if cyg_ha >= _cyg_transit_window[0] and cyg_ha <= _cyg_transit_window[1]:
            print(f"{filename} 'CygA transit' fLow={freqs[0]/1e6:.1f} MHz fHigh={freqs[-1]/1e6:.1f} MHz LST={lst:.2f} hr")
        elif lst >= _diffuse_lst_window[0] and lst <= _diffuse_lst_window[1]:
            print(f"{filename} 'diffuse' fLow={freqs[0]/1e6:.1f} MHz fHigh={freqs[-1]/1e6:.1f} MHz  LST={lst:.2f} hr")
        else:
            print(f"{filename} 'verification only' fLow={freqs[0]/1e6:.1f} MHz fHigh={freqs[-1]/1e6:.1f} MHz  LST={lst:.2f} hr")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="given an SSMIF and a collection of TBX files, figure out which are useful for calibrating the station",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
    parser.add_argument('ssmif', type=str,
                        help='station SSMIF')
    parser.add_argument('filename', type=str, nargs='+',
                        help='filename to process')
    args = parser.parse_args()
    main(args)
