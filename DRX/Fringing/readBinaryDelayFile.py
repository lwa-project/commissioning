#!/usr/bin/env python

"""
Simple script to read in a MCS binary packed DP delay file (.df) and print 
out the delays in ns.

$Rev$
$LastChangedBy$
$LastChangedDate$
"""

import os
import sys
import numpy
import struct

from lsl.common.dp import DPDtodelay
from lsl.common.stations import lwa1


def main(args):
    filename = args[0]

    # Read in the entire file
    fh = open(filename, 'rb')
    data = fh.read()
    fh.close()
    
    # Unpack the raw delays (520 unsigned short ints)
    rawDelays = struct.unpack('<520H', data)

    # Convert to delays in ns
    delays = [DPDtodelay(d) for d in rawDelays]
    
    # Report
    ants = lwa1.getAntennas()[0::2]
    
    print "Std   X [ns]    Y [ns]"
    print "----------------------"
    for i in xrange(len(ants)):
        dx, dy = delays[2*i+0], delays[2*i+1]
        print "%3i   %7.2f  %7.2f" % (ants[i].stand.id, dx, dy)


if __name__ == "__main__":
    main(sys.argv[1:])

