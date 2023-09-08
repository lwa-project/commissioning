#!/usr/bin/python3

"""
Creates delay file for one or all antennas
"""
# written by Gerald Crichton, Robert Navarro   
# March 12, 2011                                                        *
# Copyright 2011, by the California Institute of Technology. ALL RIGHTS RESERVED.

import os
import sys

zero = [0, 0]

def make_delayfile(path, ant, coarse, fine):
    delay = [coarse, fine]
    coarse = coarse.replace('.', '')
    fine = fine.replace('.', '')
    if ant.lower() == 'all':
        df_filename = 'delay_all_c%s_f%s.df' % (coarse, fine)
        dft_filename = df_filename+'t'
        filename = path + '/' + dft_filename
        file = open(path + '/' + dft_filename, 'w')
        for x in range(1, 520+1):
            file.write("%s %s\n" % (delay[0], delay[1]))
        file.close()
    else:
        df_filename = 'delay_a%s_c%s_f%s.df' % (ant, coarse, fine)
        dft_filename = df_filename+'t'
        file = open(path + '/' + dft_filename, 'w')
        ant = int(ant)
        for x in range(1, 520+1):
            if x==ant: 
                file.write("%s %s\n" % (delay[0], delay[1]))
            else:
                file.write("%s %s\n" % (zero[0], zero[1]))
        file.close()
    return [dft_filename,df_filename]
    
# take a list of delay values in ns and convert to dfile    
def list2delayfile(path, filename, dlylist):
    # convert from delay in ns to delay in fine samples
    # fine samples have a period of 1/(196000000 * 16) seconds
    finedelay_period = (1.0/(196000000*16.0)) * 1e9

    coarselist = []
    finelist = []
    fsamplelist = [0]*520
    for i in range(520):
        fsamplelist[i] = round(dlylist[i]/finedelay_period)
        coarse = fsamplelist[i] // 16
        fine   = fsamplelist[i] % 16
        if coarse > 1023:
            print("Warning: Coarse delay greater than 1023")
        coarselist.append(int(coarse))
        finelist.append(int(fine))
    df_filename = filename + '.df'
    dft_filename = df_filename+'t'
    file = open(path + '/' + dft_filename, 'w')      
    for i in range(520):
        file.write("%s %s\n" % (coarselist[i], finelist[i]))
    file.close()
    return [dft_filename,df_filename]   

########################################################
########################################################

if __name__ == '__main__':
    if len(sys.argv) < 4:
        print('usage: delay ant coarse fine')
        print('       where ant = antenna number or ''all'' ')
        exit()
        
    filename = make_delayfile('.', sys.argv[1], sys.argv[2], sys.argv[3])
    print(filename)
    
