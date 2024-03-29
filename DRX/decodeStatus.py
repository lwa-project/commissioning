#!/usr/bin/env python3

"""
Given a DP board status code, decode it into its various parts.
"""

import sys
import struct

def num2num(v):
    """
    Convert a 32-bit unsigned integer to its binary representation split 
    into blocks of four bits.  Returns the binary representation as a 
    string.
    """

    out = ""
    for i in range(32):
        if i % 4 == 0 and i != 0:
            out += " "
        out += str((v>>(31-i))&1)
    return out


def main(args):
    # Get the argument as a int
    code = int(args[0])

    # Print out the code and its binary representation
    print("Numeric Code: %i" % code)
    print("Binary Code [MSB to LSB]:  %s" % num2num(code))
    print("")

    # Decode the four status codes for each of the five FPGAs
    for i,chip in enumerate([1, 2, 3, 4, 5]):
        print("FPGA #%i" % chip)
        for j,field in enumerate(['calib', 'temp', 'vccint', 'vccaux']):
            if not (code>>(i*4+j))&1:
                status = "OK"
            else:
                status = "FAIL"
            print("%-10s %s" % (("-> %s:" % field), status))

    # Decode the overall DP power status code
    print("DP")
    if not (code>>20)&1:
        status = "OK"
    else:
        status = "FAIL"
    print("%-10s %s" % (("-> %s:" % "power"), status))


if __name__ == "__main__":
    main(sys.argv[1:])
    
