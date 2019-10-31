TBN Commissioning Scripts
=========================

checkTimetags.py
----------------
Simple script for making sure that time tags update properly in a TBN file.  The script 
not only checks time tags but also checks for dropped packets and synchronization problems.

eyeDiagram.py
-------------
Create an eye diagram for some portion of a TBN file.

eyeDiagram2.py
--------------
Look for glitches in a TBN file by fitting a quantized sine wave to the data.

tbnTimeseries.py
----------------
Script to plot out time series I/Q data from a TBN file.

tbnFileCheck.py
---------------
Script to skip through a TBN file and check for clipping on stand 10.

checkMissing.py
---------------
Look at what is coming out of the TBN ring buffer and plot of what frames are missing,
if any, as a function of time in the file.

hdfWaterfall.py
---------------
Process a TBN file into an HDF5 waterfall file where each antenna is stored as a
separate tuning.

Calibration
-----------
See the Calibration/README file.

Prototype
---------
See the Prototype/README file.

