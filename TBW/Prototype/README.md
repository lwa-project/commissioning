Prototype DP System - TBW Scipts
================================
This directory provides a variety of scripts to work with the 20 output DP
Prototype system.  This system was originally installed at LWA1 before the 
delivery of the full DP and not resides at the North Arm site.  All scripts
mentioned here uses the lwana-ssmif.txt file for mappings.

tbwSpectra.py
-------------
Basic script to plot integrated spectra from all 20 dipoles.

stationMaster.py
----------------
Script to take a TBW file and compute spectra for all dipoles, as well as the 
estimated resonance point.  This script also computes statistics about the
time domain data (average power, range, etc.).

stationMaster2.py
-----------------
Version of stationMaster.py that accepts the `ClipLevel' keyword to help clean 
up burst RFI.

smGUI.py
--------
GUI for viewing and interacting with the NPZ files created by stationMaster.py/
stationMaster2.py.

