TBF Commissioning Scripts
=========================

checkTimetags.py
----------------
Script for checking for missing frames and bad time tags in a full ADP TBF capture.  
If problems are found the error are noted and printed out.  It shouldn't be too hard
to modify this for fewer inputs.

stationMasterLite.py
--------------------
Modified version of the TBW stationMaster.py script that estimates the frequency (in 
MHz) of each dipole's resonance point and saves that information along with time-
average spectra to a NPZ file.

tbfMux.py
---------
Combine multiple single-server TBF files from the ADP triggering system into a single
file.