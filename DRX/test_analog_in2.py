#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
Set DRX beams and start a recording.
"""

# Python3 compatiability
from __future__ import print_function, division
import sys
if sys.version_info > (3,):
    xrange = range
    
import datetime
import os
import math
import socket
import string
import struct
import subprocess
import sys
try:
    import thread
except ImportError:
    import _thread as thread
import time
import gain
import delay

import getopt


def usage(exitCode=None):
    print("""test_analog_in2.py - Set DRX beams and start a recording.
    
Usage: test_analog_in2.py [OPTIONS]

Options:
-h, --help           Display this help message
-1, --freq1          Set frequency of tuning 1 in MHz (default = 12.15 MHz)
-2, --freq2          Set frequency of tuning 1 in MHz (default = 12.35 MHz)
-f, --filter         Filter code (default = 4)
-g, --gain1          DRX gain for tuning 1 (default = 5)
-a, --gain2          DRX gain for tuning 2 (default = 5)
-s, --sub-plot       Sub-slot for when the command will take place (default = 0)
-t, --time           Observation time in seconds (default = 10.000)
""")
    
    if exitCode is not None:
        sys.exit(exitCode)
    else:
        return True


def parseOptions(args):
    config = {}
    # Command line flags - default values
    config['freq1'] = 12.15e6
    config['freq2'] = 12.35e6
    config['filter'] = 4
    config['gain1'] = 5
    config['gain2'] = 5
    config['slot'] = 0
    config['time'] = 10.0
    
    # Read in and process the command line flags
    try:
        opts, args = getopt.getopt(args, "h1:2:f:g:a:s:t:", ["help", "freq1=", "freq2=", "filter=", "gain1=", "gain2=", "sub-slot=", "time="])
    except getopt.GetoptError as err:
        # Print help information and exit:
        print(str(err)) # will print something like "option -a not recognized"
        usage(exitCode=2)
    
    # Work through opts
    for opt, value in opts:
        if opt in ('-h', '--help'):
            usage(exitCode=0)
        elif opt in ('-1', '--freq1'):
            config['freq1'] = float(value)*1e6
        elif opt in ('-2', '--freq2'):
            config['freq2'] = float(value)*1e6
        elif opt in ('-f', '--filter'):
            config['filter'] = int(value)
        elif opt in ('-g', '--gain1'):
            config['gain1'] = int(value)
        elif opt in ('-a', '--gain2'):
            config['gain2'] = int(value)
        elif opt in ('-s', '--sub-slot'):
            config['slot'] = int(value)
        elif opt in ('-t', '--time'):
            config['time'] = float(value)
        else:
            assert False
    
    # Validate inputs
    assert(config['filter'] >=1 and config['filter'] <= 7)
    assert(config['freq1'] >= 10.0e6 and config['freq1'] <= 88.0e6)
    assert(config['freq2'] >= 10.0e6 and config['freq2'] <= 88.0e6)
    assert(config['slot'] >= 0 and config['slot'] <= 99)
    assert(config['time'] > 0)
    assert(config['gain1'] >= 0 and config['gain1'] <= 15 )
    assert(config['gain2'] >= 0 and config['gain2'] <= 15 )
    
    # Add in arguments
    config['args'] = args

    # Return configuration
    return config


def get_time():

    # determine current time
    dt = datetime.datetime.utcnow()
    year        = dt.year             
    month       = dt.month      
    day         = dt.day    
    hour        = dt.hour
    minute      = dt.minute
    second      = dt.second     
    millisecond = dt.microsecond / 1000

    # compute MJD         
    # adapted from http://paste.lisp.org/display/73536
    # can check result using http://www.csgnetwork.com/julianmodifdateconv.html
    a = (14 - month) // 12
    y = year + 4800 - a          
    m = month + (12 * a) - 3                    
    p = day + (((153 * m) + 2) // 5) + (365 * y)   
    q = (y // 4) - (y // 100) + (y // 400) - 32045
    mjdi = int(math.floor( (p+q) - 2400000.5))
    mjd = "%6s" % mjdi      

    # compute MPM
    mpmi = int(math.floor( (hour*3600 + minute*60 + second)*1000 + millisecond ))
    mpm = string.rjust(str(mpmi),9)

    return (mjd, mpm)


def wait_til_sec():

    # determine current time
    (mjd, mpm) = get_time()
    msecs = float(mpm)
    sec = int(msecs/1000)
    frac = msecs - (sec*1000)
    sleeptime = ((1000-frac)/1000) + 0.005
    time.sleep(sleeptime)    

    (mjd, mpm) = get_time()
    return (mjd, mpm)
###############################################################################
if __name__ == '__main__':

    config = parseOptions(sys.argv[1:])

    # Get current time
    #(mjd, mpm) = get_time()
    (mjd, mpm) = wait_til_sec()
    print('starting at mjd=%i; mpm=%i' % (int(mjd), int(mpm)))
    
    execpath = '/home/ops/JR5/src/exec' 
    schpath =  '/home/ops/JR5/src/sch'
    gfilepath = schpath + '/gfiles'
    dfilepath = schpath + '/dfiles'
    cfilepath = schpath + '/cfiles'
    
    
    # Create delays and gains files
    bgain = 1.0000
    bgain_cross = 0.0001
    gainlist = [[0,0,0,0]]*260 # initialize gain list
    for x in range(0,8):
        gainlist[x] = [bgain, bgain_cross, bgain_cross, bgain] # create list gains with no cross polarization correction
    gfiles = gain.list2gainfile(gfilepath, 'analog_gain1', gainlist)
    os.system( os.path.join(execpath,'megfg')
        + ' ' + os.path.join(gfilepath,gfiles[0])
        + ' ' + os.path.join(gfilepath,gfiles[1]))

    sdelay = 1/20.0e6/2*1e9
    delaylist = [0]*520
    for x in range(0,8):
        delaylist[x] = sdelay
    dfiles = delay.list2delayfile(dfilepath, 'analog_delay1', delaylist)
    #dfiles = delay.make_delayfile(dfilepath, 'all', '0', '0')
    os.system( os.path.join(execpath,'medfg')
        + ' ' + os.path.join(dfilepath,dfiles[0])
        + ' ' + os.path.join(dfilepath,dfiles[1]))     
    

    print('Set beam 1 delays and gains, then set DRX for beam 1')
    os.system(os.path.join(execpath,'mesix') + ' DP_ BAM "1 ' + dfiles[1] + ' ' + gfiles[1] + ' 0" today +5')
    os.system(os.path.join(execpath,'mesix') + ' DP_ BAM "1 ' + dfiles[1] + ' ' + gfiles[1] + ' 3" today +5')
    print('Set beam 1 delays and gains, then set DRX for beam 2')
    os.system(os.path.join(execpath,'mesix') + ' DP_ BAM "2 ' + dfiles[1] + ' ' + gfiles[1] + ' 1" today +5')
    os.system(os.path.join(execpath,'mesix') + ' DP_ BAM "2 ' + dfiles[1] + ' ' + gfiles[1] + ' 4" today +5')  
    print('Set beam 1 delays and gains, then set DRX for beam 3')
    os.system(os.path.join(execpath,'mesix') + ' DP_ BAM "3 ' + dfiles[1] + ' ' + gfiles[1] + ' 2" today +5')
    os.system(os.path.join(execpath,'mesix') + ' DP_ BAM "3 ' + dfiles[1] + ' ' + gfiles[1] + ' 5" today +5') 
    print('Set beam 1 delays and gains, then set DRX for beam 4')
    os.system(os.path.join(execpath,'mesix') + ' DP_ BAM "4 ' + dfiles[1] + ' ' + gfiles[1] + ' 3" today +5')
    os.system(os.path.join(execpath,'mesix') + ' DP_ BAM "4 ' + dfiles[1] + ' ' + gfiles[1] + ' 6" today +5')          

    time.sleep(1)
    # Now set up DRX 
    os.system('%s/mesix DP_ DRX "3 1 %i %i %i %i"' % (execpath, config['freq1'], 7, config['gain1'], config['slot'])) # Beam, Tuning, Freq, BW, gain, subslot
    os.system('%s/mesix DP_ DRX "3 2 %i %i %i %i"' % (execpath, config['freq2'], 6, config['gain2'], config['slot']))
    time.sleep(1)
    os.system('%s/mesix DP_ DRX "4 1 %i %i %i %i"' % (execpath, config['freq1'], 7, config['gain1'], config['slot'])) # Beam, Tuning, Freq, BW, gain, subslot
    os.system('%s/mesix DP_ DRX "4 2 %i %i %i %i"' % (execpath, config['freq2'], 6, config['gain2'], config['slot']))


    # Assume that DR is launched
    #os.system(execpath + 'mesix DR1 INI')
    #time.sleep(1)
    print('Scheduling a recording to start about 40 seconds from current time')
    (mjd, mpm) = get_time()
    new_mpm = str(int(mpm) + 40000) # Make sure both DP_MCS and DR are in time sync first. May need to do manual sync on DR
    rec_duration = str(int(round(config['time']*1000))) # time in milliseconds
    os.system(execpath + '/mesix DR1 REC "'+ mjd  + ' ' + new_mpm + ' ' + rec_duration + ' DEFAULT_DRX"')
    # now get DR1 RPT info
    os.system(execpath + '/mesix DR1 RPT DIRECTORY-COUNT')
    time.sleep(1)   
    f = subprocess.Popen(schpath + '/ms_mdre DR1 DIRECTORY-COUNT', shell=True, cwd = schpath, stdout=subprocess.PIPE).stdout
    reply = f.read().split()  
    count = reply[0] 
    print('count = '+ count)
    #os.system(execpath + '/mesix DR1 RPT DIRECTORY-ENTRY-'+ count)
    time.sleep(1)

    f = subprocess.Popen(schpath + '/ms_mdre DR1 DIRECTORY-ENTRY-'+count, shell=True, cwd = schpath, stdout=subprocess.PIPE).stdout
    reply = f.read().split()  
    fname = reply[0] 
    bytes = reply[6] # number of bytes read
    print('fname = %s; bytes read = %s' % (count,bytes)   )
    ###########################
    #time.sleep(61) # sleep a long time to make sure recording is done before continuing      
    #os.system(execpath + 'mesix DR1 RPT DIRECTORY-ENTRY-'+ count)
    #time.sleep(1)
    #f = subprocess.Popen(schpath + 'ms_mdre DR1 DIRECTORY-ENTRY-'+count, shell=True, cwd = schpath, stdout=subprocess.PIPE).stdout
    #reply = f.read().split()  
    #fname = reply[0] 
    #bytes = reply[6] # number of bytes read
    #print('fname = %s; bytes read = %s' % (count,bytes))
    ## Use temp output filename for now. Just copy immediatly
    #print('Writing a file to external disk--------------------')
    #os.system(execpath + 'mesix DR1 CPY "' + fname + ' 0 ' + bytes + ' /dev/sdg1 drx_temp.dat"')      
        
    #OK
    print('Done \n')
