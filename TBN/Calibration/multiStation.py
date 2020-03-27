"""
Temporary module that implements a parse_ssmif() function that works with
LWA-SV and LWA-1.
"""

# Python3 compatiability
from __future__ import print_function, division
import sys
if sys.version_info > (3,):
    xrange = range
    
import os
import re
import struct

from lsl.common.stations import Antenna, Stand, FEE, Cable, ARX, LWAStation, _id2name
from lsl.common import mcs, mcsADP

__version__ = "0.1"
__all__ = ['parse_ssmif',]

def __parseTextSSMIF(filename):
    """
    Given a human-readable (text) SSMIF file and return a collection of
    variables via locals() containing the files data.
    """
    
    fh = open(filename, 'r')
    
    kwdRE = re.compile(r'(?P<keyword>[A-Z_0-9]+)(\[(?P<id1>[0-9]+?)\])?(\[(?P<id2>[0-9]+?)\])?(\[(?P<id3>[0-9]+?)\])?')
    
    # Loop over the lines in the file
    for line in fh:
        line = line.replace('\n', '')
        line = line.replace('\r', '')
        if len(line) == 0 or line.isspace():
            continue
        if line[0] == '#':
            continue
        
        keywordSection, value = line.split(None, 1)
        value = value.split('#', 1)[0]
        
        mtch = kwdRE.match(keywordSection)
        keyword = mtch.group('keyword')
        
        ids = [-1, -1, -1]
        for i in xrange(3):
            try:
                ids[i] = int(mtch.group('id%i' % (i+1)))
            except TypeError:
                pass
            
        #
        # Station Data
        #
            
        if keyword == 'STATION_ID':
            idn = str(value)
            continue
        if keyword == 'GEO_N':
            lat = float(value)
            continue
        if keyword == 'GEO_E':
            lon = float(value)
            continue
        if keyword == 'GEO_EL':
            elv = float(value)
            continue
        
        #
        # Stand & Antenna Data
        #
        
        if keyword == 'N_STD':
            nStand = int(value)
            
            stdPos = [[0.0, 0.0, 0.0] for n in xrange(nStand)]
            stdAnt = [n/2+1 for n in xrange(2*nStand)]
            stdOrie = [n % 2 for n in xrange(2*nStand)]
            stdStat = [3 for n in xrange(2*nStand)]
            stdTheta = [0.0 for n in xrange(2*nStand)]
            stdPhi = [0.0 for n in xrange(2*nStand)]
            
            stdDesi = [1 for x in xrange(2*nStand)]
            
            continue
            
        if keyword == 'STD_LX':
            stdPos[ids[0]-1][0] = float(value)
            continue
        if keyword == 'STD_LY':
            stdPos[ids[0]-1][1] = float(value)
            continue
        if keyword == 'STD_LZ':
            stdPos[ids[0]-1][2] = float(value)
            continue
        
        if keyword == 'ANT_STD':
            stdAnt[ids[0]-1] = int(value)
            continue
        
        if keyword == 'ANT_ORIE':
            stdOrie[ids[0]-1] = int(value)
            continue
        
        if keyword == 'ANT_STAT':
            stdStat[ids[0]-1] = int(value)
            continue
        
        if keyword == 'ANT_THETA':
            stdTheta[ids[0]-1] = float(value)
            continue
        
        if keyword == 'ANT_PHI':
            stdPhi[ids[0]-1] = float(value)
            continue
        
        #
        # FEE, Cable, & SEP Data
        #
        
        if keyword == 'N_FEE':
            nFee = int(value)
            
            feeID = ["UNK" for n in xrange(nFee)]
            feeStat = [3 for n in xrange(nFee)]
            feeDesi = [1 for n in xrange(nFee)]
            feeGai1 = [35.7 for n in xrange(nFee)]
            feeGai2 = [35.7 for n in xrange(nFee)]
            feeAnt1 = [2*n+1 for n in xrange(nFee)]
            feeAnt2 = [2*n+2 for n in xrange(nFee)]
            
            continue
            
        if keyword == 'FEE_ID':
            feeID[ids[0]-1] = value
            continue
        
        if keyword == 'FEE_STAT':
            feeStat[ids[0]-1] = int(value)
            continue
            
        if keyword == 'FEE_DESI':
            feeDesi[ids[0]-1] = int(value)
            continue
            
        if keyword == 'FEE_GAI1':
            feeGai1[ids[0]-1] = float(value)
            continue
        if keyword == 'FEE_GAI2':
            feeGai2[ids[0]-1] = float(value)
            continue
            
        if keyword == 'FEE_ANT1':
            feeAnt1[ids[0]-1] = int(value)
            continue
        if keyword == 'FEE_ANT2':
            feeAnt2[ids[0]-1] = int(value)
            continue
        
        
        if keyword == 'N_RPD':
            nRPD = int(value)
            
            rpdID = ['UNK' for n in xrange(nRPD)]
            rpdStat = [3 for n in xrange(nRPD)]
            rpdLeng = [0.0 for n in xrange(nRPD)]
            rpdVF = [83.0 for n in xrange(nRPD)]
            rpdDD = [2.4 for n in xrange(nRPD)]
            rpdA0 = [0.00428 for n in xrange(nRPD)]
            rpdA1 = [0.00000 for n in xrange(nRPD)]
            rpdFre = [10e6 for n in xrange(nRPD)]
            rpdStr = [1.0 for n in xrange(nRPD)]
            rpdDesi = [1 for n in xrange(nRPD)]
            rpdAnt = [n+1 for n in xrange(nRPD)]
            
            continue
        
        if keyword == 'RPD_ID':
            rpdID[ids[0]-1] = value
            continue
        
        if keyword == 'RPD_STAT':
            rpdStat[ids[0]-1] = int(value)
            continue
        
        if keyword == 'RPD_LENG':
            rpdLeng[ids[0]-1] = float(value)
            continue
        
        if keyword == 'RPD_VF':
            if ids[0] == -1:
                rpdVF = [float(value) for n in xrange(nRPD)]
            else:
                rpdVF[ids[0]-1] = float(value)
            continue
        if keyword == 'RPD_DD':
            if ids[0] == -1:
                rpdDD = [float(value) for n in xrange(nRPD)]
            else:
                rpdDD[ids[0]-1] = float(value)
            continue
        
        if keyword == 'RPD_A0':
            if ids[0] == -1:
                rpdA0 = [float(value) for n in xrange(nRPD)]
            else:
                rpdA0[ids[0]-1] = float(value)
            continue
        
        if keyword == 'RPD_A1':
            if ids[0] == -1:
                rpdA1 = [float(value) for n in xrange(nRPD)]
            else:
                rpdA1[ids[0]-1] = float(value)
            continue
        
        if keyword == 'RPD_FREF':
            if ids[0] == -1:
                rpdFre = [float(value) for n in xrange(nRPD)]
            else:
                rpdFre[ids[0]-1] = float(value)
            continue
        
        if keyword == 'RPD_STR':
            if ids[0] == -1:
                rpdStr = [float(value) for n in xrange(nRPD)]
            else:
                rpdStr[ids[0]-1] = float(value)
            continue
        
        if keyword == 'RPD_DESI':
            rpdDesi[ids[0]-1] = value
            continue
        
        if keyword == 'RPD_ANT':
            rpdAnt[ids[0]-1] = int(value)
            continue
        
        
        if keyword == 'N_SEP':
            nSEP = int(value)
            
            sepCbl = ['UNK' for n in xrange(nSEP)]
            sepLeng = [0.0 for n in xrange(nSEP)]
            sepDesi = [1 for n in xrange(nSEP)]
            sepGain = [0.0 for n in xrange(nSEP)]
            sepAnt = [n+1 for n in xrange(nSEP)]
            
            continue
        
        if keyword == 'SEP_CABL':
            sepCbl[ids[0]-1] = value
            continue
        
        if keyword == 'SEP_LENG':
            sepLeng[ids[0]-1] = float(value)
            continue
        
        if keyword == 'SEP_DESI':
            sepDesi[ids[0]-1] = int(value)
            continue
        
        if keyword == 'SEP_GAIN':
            sepGain[ids[0]-1] = float(value)
            continue
        
        if keyword == 'SEP_ANT':
            sepAnt[ids[0]-1] = int(value)
            continue
        
        #
        # ARX (ARB) Data
        #
        
        if keyword == 'N_ARB':
            nARX = int(value)
            
            arxID = ["UNK" for n in xrange(nARX)]
            arxSlot = [0 for n in xrange(nARX)]
            arxDesi = [0 for n in xrange(nARX)]
            arxRack = [0 for n in xrange(nARX)]
            arxPort = [0 for n in xrange(nARX)]
            
            continue
        
        if keyword == 'N_ARBCH':
            nchanARX = int(value)
            
            arxStat = [[3 for c in xrange(nchanARX)] for n in xrange(nARX)]
            arxAnt = [[n*nchanARX+c+1 for c in xrange(nchanARX)] for n in xrange(nARX)]
            arxIn = [["UNK" for c in xrange(nchanARX)] for n in xrange(nARX)]
            arxOut = [["UNK" for c in xrange(nchanARX)] for n in xrange(nARX)]
            
            continue
        
        if keyword == 'ARB_ID':
            arxID[ids[0]-1] = value
            continue
        
        if keyword == 'ARB_SLOT':
            try:
                arxSlot[ids[0]-1] = int(value)
            except ValueError:
                arxSlot[ids[0]-1] = value
        
        if keyword == 'ARB_DESI':
            arxDesi[ids[0]-1] = int(value)
            continue
        
        if keyword == 'ARB_RACK':
            arxRack[ids[0]-1] = int(value)
            continue
        
        if keyword == 'ARB_PORT':
            arxRack[ids[0]-1] = int(value)
            continue
        
        if keyword == 'ARB_STAT':
            arxStat[ids[0]-1][ids[1]-1] = int(value)
            continue
        
        if keyword == 'ARB_ANT':
            arxAnt[ids[0]-1][ids[1]-1] = int(value)
            continue
        
        if keyword == 'ARB_IN':
            arxIn[ids[0]-1][ids[1]-1] = value
            continue
        
        if keyword == 'ARB_OUT':
            arxOut[ids[0]-1][ids[1]-1] = value
            continue
        
        #
        # DP 1 & 2 Data - LWA1
        #
        
        if keyword == 'N_DP1':
            nDP1 = int(value)
            
            dp1ID = ["UNK" for n in xrange(nDP1)]
            dp1Slot = [0 for n in xrange(nDP1)]
            dp1Desi = [1 for n in xrange(nDP1)]
            
            continue
        
        if keyword == 'N_DP1CH':
            nchanDP1 = int(value)
            
            dp1Stat = [[3 for c in xrange(nchanDP1)] for n in xrange(nDP1)]
            dp1InR = [["UNK" for c in xrange(nchanDP1)] for n in xrange(nDP1)]
            dp1InC = [["UNK" for c in xrange(nchanDP1)] for n in xrange(nDP1)]
            dp1Ant = [[n*nchanDP1+c+1 for c in xrange(nchanDP1)] for n in xrange(nDP1)]
            
            continue
        
        if keyword == 'DP1_ID':
            dp1ID[ids[0]-1] = value
            continue
        
        if keyword == 'DP1_SLOT':
            dp1Slot[ids[0]-1] = value
            continue
        
        if keyword == 'DP1_DESI':
            dp1Desi[ids[0]-1] = int(value)
            continue
        
        if keyword == 'DP1_STAT':
            dp1Stat[ids[0]-1][ids[1]-1] = int(value)
            continue
        
        if keyword == 'DP1_INR':
            dp1InR[ids[0]-1][ids[1]-1] = value
            continue
        
        if keyword == 'DP1_INC':
            dp1InC[ids[0]-1][ids[1]-1] = value
            continue
        
        if keyword == 'DP1_ANT':
            dp1Ant[ids[0]-1][ids[1]-1] = int(value)
        
        
        if keyword == 'N_DP2':
            nDP2 = int(value)
            
            dp2ID = ["UNK" for n in xrange(nDP2)]
            dp2Slot = ["UNK" for n in xrange(nDP2)]
            dp2Stat = [3 for n in xrange(nDP2)]
            dp2Desi = [1 for n in xrange(nDP2)]
            
            continue
        
        if keyword == 'DP2_ID':
            dp2ID[ids[0]-1] = value
            continue
        
        if keyword == 'DP2_SLOT':
            dp2Slot[ids[0]-1] = value
            continue
        
        if keyword == 'DP2_STAT':
            dp2Stat[ids[0]-1] = int(value)
            continue
        
        if keyword == 'DP2_DESI':
            dp2Desi[ids[0]-1] = int(value)
            continue
        
        #
        # ROACH & Server Data - LWA-SV
        #
        
        if keyword == 'N_ROACH':
            nRoach = int(value)
            
            roachID = ["UNK" for n in xrange(nRoach)]
            roachSlot = [0 for n in xrange(nRoach)]
            roachDesi = [1 for n in xrange(nRoach)]
            
            continue
        
        if keyword == 'N_ROACHCH':
            nchanRoach = int(value)
            
            roachStat = [[3 for c in xrange(nchanRoach)] for n in xrange(nRoach)]
            roachInR = [["UNK" for c in xrange(nchanRoach)] for n in xrange(nRoach)]
            roachInC = [["UNK" for c in xrange(nchanRoach)] for n in xrange(nRoach)]
            roachAnt = [[n*nchanRoach+c+1 for c in xrange(nchanRoach)] for n in xrange(nRoach)]
            
            continue
        
        if keyword == 'ROACH_ID':
            roachID[ids[0]-1] = value
            continue
        
        if keyword == 'ROACH_SLOT':
            roachSlot[ids[0]-1] = value
            continue
        
        if keyword == 'ROACH_DESI':
            roachDesi[ids[0]-1] = int(value)
            continue
        
        if keyword == 'ROACH_STAT':
            roachStat[ids[0]-1][ids[1]-1] = int(value)
            continue
        
        if keyword == 'ROACH_INR':
            roachInR[ids[0]-1][ids[1]-1] = value
            continue
        
        if keyword == 'ROACH_INC':
            roachInC[ids[0]-1][ids[1]-1] = value
            continue
        
        if keyword == 'ROACH_ANT':
            roachAnt[ids[0]-1][ids[1]-1] = int(value)
        
        
        if keyword == 'N_SERVER':
            nServer = int(value)
            
            serverID = ["UNK" for n in xrange(nServer)]
            serverSlot = ["UNK" for n in xrange(nServer)]
            serverStat = [3 for n in xrange(nServer)]
            serverDesi = [1 for n in xrange(nServer)]
            
            continue
        
        if keyword == 'SERVER_ID':
            serverID[ids[0]-1] = value
            continue
        
        if keyword == 'SERVER_SLOT':
            serverSlot[ids[0]-1] = value
            continue
        
        if keyword == 'SERVER_STAT':
            serverStat[ids[0]-1] = int(value)
            continue
        
        if keyword == 'SERVER_DESI':
            serverDesi[ids[0]-1] = int(value)
            continue
        
        #
        # DR Data
        #
        
        if keyword == 'N_DR':
            nDR = int(value)
            
            drStat = [0 for n in xrange(nDR)]
            drID = ["UNK" for n in xrange(nDR)]
            drShlf = [0 for n in xrange(nDR)]
            drPC = ["UNK" for n in xrange(nDR)]
            drDP = [0 for n in xrange(nDR)]
            
            continue
        
        if keyword == 'DR_STAT':
            drStat[ids[0]-1] = int(value)
            continue
        
        if keyword == 'DR_ID':
            drID[ids[0]-1] = value
            continue
        
        if keyword == 'DR_SHLF':
            drShlf[ids[0]-1] = int(value)
            continue
        
        if keyword == 'DR_PC':
            drPC[ids[0]-1] = value
            continue
        
        if keyword == 'DR_DP':
            drDP[ids[0]-1] = int(value)
            continue
        
    fh.close()
    
    return locals()


def __parseBinarySSMIF(filename):
    """
    Given a binary packed SSMIF file and return a collection of
    variables via locals() containing the files data.
    """
    
    fh = open(filename, 'rb')
    
    # Read in the first four bytes to get the version code and go from there
    version = fh.read(4)
    version = struct.unpack('<i', version)[0]
    fh.seek(0)
    
    if version in (8,):
        ## ADP
        mode = mcsADP
    else:
        ## DP
        mode = mcs
    bssmif = mode.parse_c_struct(mode.SSMIF_STRUCT, char_mode='int', endianness='little')
    bsettings = mode.parse_c_struct(mode.STATION_SETTINGS_STRUCT, endianness='little')
    
    fh.readinto(bssmif)
    
    #
    # Station Data
    #
    idn = [chr(i) for i in bssmif.sStationID]
    idn = ''.join([i for i in idn if i != '\x00'])
    lat = bssmif.fGeoN
    lon = bssmif.fGeoE
    elv = bssmif.fGeoEl
    
    #
    # Stand & Antenna Data
    #
    stdPos   = [list(i) for i in zip(bssmif.fStdLx, bssmif.fStdLy, bssmif.fStdLz)]
    stdAnt   = list(bssmif.iAntStd)
    stdOrie  = list(bssmif.iAntOrie)
    stdStat  = list(bssmif.iAntStat)
    stdTheta = list(bssmif.fAntTheta)
    stdPhi   = list(bssmif.fAntPhi)
    stdDesi  = list(bssmif.eAntDesi)
    
    #
    # FEE, Cable, & SEP Data
    #
    feeID   = flat_to_multi([chr(i) for i in bssmif.sFEEID], *bssmif.dims['sFEEID'])
    feeID   = [''.join([k for k in i if k != '\x00']) for i in feeID]
    feeStat = list(bssmif.iFEEStat)
    feeDesi = list(bssmif.eFEEDesi)
    feeGai1 = list(bssmif.fFEEGai1)
    feeGai2 = list(bssmif.fFEEGai2)
    feeAnt1 = list(bssmif.iFEEAnt1)
    feeAnt2 = list(bssmif.iFEEAnt2)
    
    rpdID   = flat_to_multi([chr(i) for i in bssmif.sRPDID], *bssmif.dims['sRPDID'])
    rpdID   = [''.join([k for k in i if k != '\x00']) for i in rpdID]
    rpdStat = list(bssmif.iRPDStat)
    rpdDesi = list(bssmif.eRPDDesi)
    rpdLeng = list(bssmif.fRPDLeng)
    rpdVF   = list(bssmif.fRPDVF)
    rpdDD   = list(bssmif.fRPDDD)
    rpdA0   = list(bssmif.fRPDA0)
    rpdA1   = list(bssmif.fRPDA1)
    rpdFre  = list(bssmif.fRPDFref)
    rpdStr  = list(bssmif.fRPDStr)
    rpdAnt  = list(bssmif.iRPDAnt)
    
    sepCbl  = flat_to_multi([chr(i) for i in bssmif.sSEPCabl], *bssmif.dims['sSEPCabl'])
    sepCbl  = [''.join([k for k in i if k != '\x00']) for i in sepCbl]
    sepLeng = list(bssmif.fSEPLeng)
    sepDesi = list(bssmif.eSEPDesi)
    sepGain = list(bssmif.fSEPGain)
    sepAnt  = list(bssmif.iSEPAnt)
    
    #
    # ARX (ARB) Data
    #
    nchanARX = bssmif.nARBCH
    arxID    = flat_to_multi([chr(i) for i in bssmif.sARBID], *bssmif.dims['sARBID'])
    arxID    = [''.join([k for k in i if k != '\x00']) for i in arxID]
    arxSlot  = list(bssmif.iARBSlot)
    arxDesi  = list(bssmif.eARBDesi)
    arxRack  = list(bssmif.iARBRack)
    arxPort  = list(bssmif.iARBPort)
    arxStat  = flat_to_multi(bssmif.eARBStat, *bssmif.dims['eARBStat'])
    arxAnt   = flat_to_multi(bssmif.iARBAnt, *bssmif.dims['iARBAnt'])
    arxIn    = flat_to_multi([chr(i) for i in bssmif.sARBIN], *bssmif.dims['sARBIN'])
    arxIn    = [[''.join(i) for i in j] for j in arxIn]
    arxOut   = flat_to_multi([chr(i) for i in bssmif.sARBOUT], *bssmif.dims['sARBOUT'])
    arxOUt   = [[''.join(i) for i in j] for j in arxOut]
    
    try:
        #
        # DP 1 & 2 Data
        #
        dp1ID   = flat_to_multi([chr(i) for i in bssmif.sDP1ID], *bssmif.dims['sDP1ID'])
        dp1ID   = [''.join([k for k in i if k != '\x00']) for i in dp1ID]
        dp1Slot = flat_to_multi([chr(i) for i in bssmif.sDP1Slot], *bssmif.dims['sDP1Slot'])
        dp1Slot = [''.join([k for k in i if k != '\x00']) for i in dp1Slot]
        dp1Desi = list(bssmif.eDP1Desi)
        dp1Stat = list(bssmif.eDP1Stat)
        dp1InR  = flat_to_multi([chr(i) for i in bssmif.sDP1INR], *bssmif.dims['sDP1INR'])
        dp1InR  = [[''.join([k for k in i if k != '\x00']) for i in j] for j in dp1InR]
        dp1InC  = flat_to_multi([chr(i) for i in bssmif.sDP1INC], *bssmif.dims['sDP1INC'])
        dp1InC  = [[''.join([k for k in i if k != '\x00']) for i in j] for j in dp1InC]
        dp1Ant  = flat_to_multi(bssmif.iDP1Ant, *bssmif.dims['iDP1Ant'])
        
        dp2ID   = flat_to_multi([chr(i) for i in bssmif.sDP2ID], *bssmif.dims['sDP2ID'])
        dp2ID   = [''.join([k for k in i if k != '\x00']) for i in dp2ID]
        dp2Slot = flat_to_multi([chr(i) for i in bssmif.sDP2Slot], *bssmif.dims['sDP2Slot'])
        dp2Slot = [''.join([k for k in i if k != '\x00']) for i in dp2Slot]
        dp2Stat = list(bssmif.eDP2Stat)
        dp2Desi = list(bssmif.eDP2Desi)
    except AttributeError:
        #
        # ROACH & Server Data
        #
        roachID   = flat_to_multi([chr(i) for i in bssmif.sRoachID], *bssmif.dims['sRoachID'])
        roachID   = [''.join([k for k in i if k != '\x00']) for i in roachID]
        roachSlot = flat_to_multi([chr(i) for i in bssmif.sRoachSlot], *bssmif.dims['sRoachSlot'])
        roachSlot = [''.join([k for k in i if k != '\x00']) for i in roachSlot]
        roachDesi = list(bssmif.eRoachDesi)
        roachStat = list(bssmif.eROoachStat)
        roachInR  = flat_to_multi([chr(i) for i in bssmif.sRoachINR], *bssmif.dims['sRoachINR'])
        roachInR  = [[''.join([k for k in i if k != '\x00']) for i in j] for j in roachInR]
        roachInC  = flat_to_multi([chr(i) for i in bssmif.sRoachINC], *bssmif.dims['sRoachINC'])
        roachInC  = [[''.join([k for k in i if k != '\x00']) for i in j] for j in roachInC]
        roachAnt  = flat_to_multi(bssmif.iRoachAnt, *bssmif.dims['iRoachAnt'])
        
        serverID   = flat_to_multi([chr(i) for i in bssmif.sServerID], *bssmif.dims['sServerID'])
        serverID   = [''.join([k for k in i if k != '\x00']) for i in serverID]
        serverSlot = flat_to_multi([chr(i) for i in bssmif.sServerSlot], *bssmif.dims['sServerSlot'])
        serverSlot = [''.join([k for k in i if k != '\x00']) for i in serverSlot]
        serverStat = list(bssmif.eServerStat)
        serverDesi = list(bssmif.eServerDesi)
        
    #
    # DR Data
    #
    drStat = list(bssmif.eDRStat)
    drID   = flat_to_multi([chr(i) for i in bssmif.sDRID], *bssmif.dims['sDRID'])
    drID   = [''.join([k for k in i if k != '\x00']) for i in drID]
    drShlf = [0 for i in xrange(bssmif.nDR)]
    drPC   = flat_to_multi([chr(i) for i in bssmif.sDRPC], *bssmif.dims['sDRPC'])
    drPC   = [''.join([k for k in i if k != '\x00']) for i in drPC]
    drDP   = list(bssmif.iDRDP)
    
    fh.readinto(bsettings)
    
    fh.close()
    
    return locals()


def parse_ssmif(filename):
    """
    Given a SSMIF file, return a fully-filled LWAStation instance.  This function
    supports both human-readable files (filenames with '.txt' extensions) or 
    binary packed files (filenames with '.dat' extensions).
    """
    
    # Find out if we have a .txt or .dat file and process accordingly
    base, ext = os.path.splitext(filename)
    
    # Read in the ssmif to a dictionary of variables
    if ext == '.dat':
        ssmifDataDict = __parseBinarySSMIF(filename)
    elif ext == '.txt':
        ssmifDataDict = __parseTextSSMIF(filename)
    else:
        raise ValueError("Unknown file extension '%s', cannot tell if it is text or binary" % ext)
    
    # Unpack the dictionary into the current variable scope
    for k in ssmifDataDict.keys():
        exec("%s = ssmifDataDict['%s']" % (k, k))
    
    # Build up a list of Stand instances and load them with data
    i = 1
    stands = []
    for pos in stdPos:
        stands.append(Stand(i, *pos))
        i += 1
    
    # Build up a list of FEE instances and load them with data
    i = 1
    fees = []
    for id,gain1,gain2,stat in zip(feeID, feeGai1, feeGai2, feeStat):
        fees.append(FEE(id, i, gain1=gain1, gain2=gain2, status=stat))
        i += 1
    
    # Build up a list of Cable instances and load them with data
    i = 1
    cables = []
    for id,length,vf,dd,a0,a1,stretch,aFreq in zip(rpdID, rpdLeng, rpdVF, rpdDD, rpdA0, rpdA1, rpdStr, rpdFre):
        cables.append(Cable(id, length, vf=vf/100.0, dd=float(dd)*1e-9, a0=a0, a1=a1, aFreq=aFreq, stretch=stretch))
        i += 1
    
    # Build up a list of Antenna instances and load them with antenna-level
    # data
    i = 1
    antennas = []
    for ant,pol,theta,phi,stat in zip(stdAnt, stdOrie, stdTheta, stdPhi, stdStat):
        antennas.append(Antenna(i, stand=stands[ant-1], pol=pol, theta=theta, phi=phi, status=stat))
        i += 1
    
    # Associate FEEs with Antennas and set the FEE port numbers
    i = 1
    for fee,ant in zip(fees, feeAnt1):
        antennas[ant-1].fee = fee
        antennas[ant-1].feePort = 1
        i += 1
    i = 1
    for fee,ant in zip(fees, feeAnt2):
        antennas[ant-1].fee = fee
        antennas[ant-1].feePort = 2
        i += 1
    
    # Associate Cables with Antennas
    i = 1
    for cbl,ant in zip(cables, rpdAnt):
        antennas[ant-1].cable = cbl
        i += 1
    
    # Associate ARX boards/channels with Antennas
    for i in xrange(len(arxAnt)):
        for j in xrange(len(arxAnt[i])):
            ant = arxAnt[i][j]
            if ant == 0 or ant > 520:
                continue
            
            boardID = arxID[i]
            channel = j + 1
            antennas[ant-1].arx = ARX(boardID, channel=channel, aspChannel=i*nchanARX + j + 1, input=arxIn[i][j], output=arxOut[i][j])
            
    try:
        # Associate DP 1 board and digitizer numbers with Antennas - DP1 boards are 2-14 and 16-28 
        # with DP2 boards at 1 and 15.
        i = 1
        j = 1
        for brd,inp in zip(dp1Ant,dp1InR):
            for ant,con in zip(brd,inp):
                antennas[ant-1].board = i + 1 + (i/14)
                antennas[ant-1].digitizer = j
                antennas[ant-1].input = con
                j += 1
            i += 1
    except NameError:
        # Associate ROACH board and digitizer numbers with Antennas.
        i = 1
        j = 1
        for brd,inp in zip(roachAnt,roachInR):
            for ant,con in zip(brd,inp):
                antennas[ant-1].board = i
                antennas[ant-1].digitizer = j
                antennas[ant-1].input = con
                j += 1
            i += 1
            
    # Build a Station
    try:
        station = LWAStation(_id2name[idn], lat, lon, elv, id=idn, antennas=antennas)
    except KeyError:
        station = LWAStation('New LWA Station', lat, lon, elv, id=idn, antennas=antennas)
    
    # And return it
    return station
