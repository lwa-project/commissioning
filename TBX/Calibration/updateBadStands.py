#!/usr/bin/env python3

import os
import sys
import argparse
from datetime import datetime, timezone

from lsl.common.stations import parse_ssmif
from lsl.logger import LSL_LOGGER

from typing import List, Tuple


def parse_stand_string(stand_csv: str) -> List[Tuple[int,int]]:
    """
    Given a comma seprated lists of stands/pol, convert that list into a
    collection of two-element tuples of (stand, pol).  The decode used is:
     * "123X" -> (123, 0)
     * "256Y" -> (256, 1)
     * "96" -> (96, 0) and (96, 1)
    """
    
    if stand_csv is None:
        return []
        
    output = []
    for s in stand_csv.upper().split(','):
        std = int(s.replace('X', '').replace('Y', ''), 10)
        pols = []
        for i,p in enumerate(['X', 'Y']):
            if s.find(p) != -1:
                pols.append(i)
        if not pols:
            pols = [0, 1]
            
        for p in pols:
            output.append((std,p))
            
    return output


def main(args):
    #
    # Load in the data
    #
    ssmifContents = open(args.filename, 'r').readlines()
    site     = parse_ssmif(args.filename)
    
    #
    # Load in the current status
    #
    ant_stat = [3,]*len(site.antennas)
    ant_comm = ['',]*len(site.antennas)
    for line in ssmifContents:
        if line.startswith('ANT_STAT'):
            start = line.find('[')
            stop  = line.find(']')
            try:
                value, toSave = line.split('#', 1)
                toSave = toSave.strip()
            except ValueError:
                value = line
                toSave = ''
            
            antID = int(line[start+1:stop], 10)
            curStatus = int(value.split(']')[1].strip(), 10)
            a = site.antennas[antID-1]
            
            LSL_LOGGER.debug(f"Found exiting status flag of '{curStatus}' for {a.stand.id}{'Y' if a.pol else 'X'}")
            ant_stat[antID-1] = curStatus
            ant_comm[antID-1] = toSave
            
    #
    # Update the antenna status flags
    #
    good = parse_stand_string(args.good)
    sus = parse_stand_string(args.suspect)
    bad = parse_stand_string(args.bad)
    
    for a in site.antennas:
        sp = (a.stand.id, a.pol)
        if sp in good:
            a.status = 3
        elif sp in sus:
            a.status = 2
        elif sp in bad:
            a.status = 1
            
        if a.status != ant_stat[a.id-1]:
            LSL_LOGGER.info(f"Status flag changes for '{ant_stat[a.id-1]}' to '{a.status}' for {a.stand.id}{'Y' if a.pol else 'X'}")
            ant_stat[a.id-1] = a.status
            ant_comm[a.id-1] = args.comment
            
    #
    # Final results
    #
    if args.output is not None:
        LSL_LOGGER.info(f"Writing updated SSMIF to '{args.output}'")
        fh = open(args.output, 'w')
    else:
        fh = sys.stdout
        
    for line in ssmifContents:
        if line.startswith('# *** ANT_STAT[antenna_id] goes here:'):
            fh.write(line)
            for i,(s,c) in enumerate(zip(ant_stat, ant_comm)):
                if s == 3:
                    continue
                    
                pad = ''
                if i+1 < 10:
                    pad = '  '
                elif i+1 < 100:
                    pad = ' '
                if c != '':
                    line = f"ANT_STAT[{i+1}]{pad} {s} # {c}\n"
                else:
                    line = f"ANT_STAT[{i+1}]{pad} {s}\n"
                fh.write(line)
        elif line.startswith('ANT_STAT'):
            continue
        else:
            fh.write(line)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Update an SSMIF to refresh the list of good, suspect, and bad stands',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
    parser.add_argument('filename', type=str,
                        help='SSMIF to update')
    parser.add_argument('-g', '--good', type=str,
                        help='command separated list of stands/pols to mark as good')
    parser.add_argument('-s', '--suspect', type=str,
                        help='command separated list of stands/pols to mark as suspect')
    parser.add_argument('-b', '--bad', type=str,
                        help='command separated list of stands/pols to mark as bad')
    parser.add_argument('-c', '--comment', type=str, default=f"Status updated on {datetime.now(tz=timezone.utc).strftime('%Y/%m/%d')}",
                        help='comment to tag updated ANT_STAT entries with')
    parser.add_argument('-o', '--output', type=str,
                        help='write output to the specified filename instead of the screen')
    args = parser.parse_args()
    main(args)
