#!/usr/bin/python3

"""
Creates gain file for one or all antennas"""

# written by Gerald Crichton, Robert Navarro   
# March 12, 2011                                                        *
# Copyright 2011, by the California Institute of Technology. ALL RIGHTS RESERVED.  

import os
import sys

zero = [0, 0, 0, 0]

def make_gainfile(path,stand, xx, xy, yx, yy):
    gain = [xx, xy, yx, yy]
    xx = xx.replace('.', '')
    xy = xy.replace('.', '')
    yx = yx.replace('.', '')
    yy = yy.replace('.', '')
    if stand.lower() == 'all':
        gf_filename = 'gain_all_%s_%s_%s_%s.gf' % (xx, xy, yx, yy)
        gft_filename = gf_filename + 't'
        file = open(path+'/'+gft_filename, 'w')
        for x in range(1, 260+1):
            file.write("%s %s %s %s\n" % (gain[0], gain[1], gain[2], gain[3]))
        file.close()
    else:
        gf_filename = 'gain_s%s_%s_%s_%s_%s.gf' % (stand, xx, xy, yx, yy)
        gft_filename = gf_filename + 't'
        file = open(path+'/'+gft_filename, 'w')
        stand = int(stand)
        for x in range(1, 260+1):
            if x==stand: 
                file.write("%s %s %s %s\n" % (gain[0], gain[1], gain[2], gain[3]))
            else:
                file.write("%s %s %s %s" % (zero[0], zero[1], zero[2], zero[3]))
            file.close()
    return [gft_filename,gf_filename]

#take a list of gain values  and convert to gfile    
def list2gainfile(path, filename, gainlist):
    gainlist_len = len(gainlist)
    if gainlist_len < 260:
        print('Gain list does not cover all 260 stands')
        # build list up to 260
        for i in range(260 - gainlist_len):
            gainlist.append([0,0,0,0])    
    for i in range(260):
        if len(gainlist[i]) != 4:
            print('stand %d does not have all four gain values')
            gainlist[i] = [0,0,0,0]
    gf_filename = filename + '.gf'
    gft_filename = gf_filename+'t'
    file = open(path + '/' + gft_filename, 'w')      
    for i in range(260):
        file.write("%s %s %s %s\n" % (gainlist[i][0], gainlist[i][1], gainlist[i][2], gainlist[i][3]))
    file.close()
    return [gft_filename,gf_filename]   
    return  1          

########################################################
########################################################

if __name__ == '__main__':
    if len(sys.argv) < 6:
        print('usage: gain stand xx xy yx yy')
        print('       where stand = stand number or ''all'' ')
        exit()
        
    filename = make_gainfile('.', sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
    print(filename)
    
