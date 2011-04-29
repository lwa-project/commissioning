#!/bin/bash

#
# checkDRSU.sh
# 
# Description: 
# Shell script to run on Linux machines with attached and powered-on DRSU to 
# report the firmware version of Seagate drives with model number "ST31000525SV".  
# This is useful for identifying disks with the newer (and less reliable) firmware 
# version "CV12".
#
# Usage:
# sudo ./checkDRSU.sh
#
# Outputs:
# List of firmware versions for all attached drives with the above model number.
#

cat /proc/scsi/scsi | grep ST31000525SV | sed -e 's/.*Rev:\s*//g;'

