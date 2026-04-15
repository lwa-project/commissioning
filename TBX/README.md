TBX Commissioning Scripts
=========================

checkTimetags.py
----------------
Script for checking for missing frames and bad time tags in a full NDP TBX capture.
If problems are found the error are noted and printed out.  It shouldn't be too hard
to modify this for fewer inputs.

stationMasterLite.py
--------------------
Modified version of the old TBW stationMaster.py script that estimates the frequency
(in MHz) of each dipole's resonance point and saves that information along with time-
average spectra to a NPZ file.

smGUI.py
--------
GUI that interfaces with a NPZ file created by stationMasterLite.py that makes looking
for problems with dipoles/mappings/etc. a point-and-click exercise.

tbxMux.py
---------
Combine multiple single-server TBX files from the NDP triggering system into a single
file.

Calibration
-----------
Collection of station calibration scripts.
