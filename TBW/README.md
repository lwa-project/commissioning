TBW Commissioning Scripts
=========================

checkTimetags.py
----------------
Script for checking for missing frames and bad time tags in a full DP TBW capture.  
If problems are found the error are noted and printed out.  It shouldn't be too hard
to modify this for fewer inputs.

tbw2hdf.py
----------
Export select stands from a TBW file to HDF5.

stationMaster.py
----------------
Modified version of tbwSpectra.py that is included with LSL (version >= 0.4.0) that 
estimated the frequency (in MHz) of each dipole's resonance point and saves that 
information along with time-average spectra to a NPZ file.

tinyMaster.py
-------------
Modified version of stationMaster.py that uses Numpy memory maps to work on systems with
less than 16 GB of memory.

stationMaster2.py
-----------------
Modified version of stationMaster.py that uses the 'ClipLevel' keyword in LSL >= 0.4.2
to blank impulsive RFI events.

tinyMaster2.py
--------------
Modified version of stationMaster2.py that uses Numpy memory maps to work on systems with
less than 16 GB of memory.

healthCheck2stationMaster.py
----------------------------
Convert a TBW health check file into a .npz file that works with smGUI.py.

smGUI.py
--------
GUI that interfaces with a NPZ file created by stationMaster.py that makes looking for
problems with dipoles/mappings/etc. a point-and-click exercise.

rfiCheck.py
-----------
Script to take a single TBW capture and create a RFI-centered HDF5 file for stands 1, 10, 54, 
248, 251, and 258 (the outlier).  These stands correspond to the four corners of the array, the
center, and the outlier.  The HDF5 contains values for the spectral kurtosis estimated from
the data and various statistics about the timeseries (mean, std. dev., percentiles, etc.)

plotRFI.py
----------
Script to take the output of rfiCheck.py and make a series of plots to look at.


