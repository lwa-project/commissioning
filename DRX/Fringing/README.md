DRX Commissioning Scripts - Fringing
====================================

fringeSDF.py
------------
Read in SSMIF file and create a set a SDF that puts a single dipole or beam on the 
X pol and the outlier on the other.

fringeDipole.py
---------------
Script to fringe DRX files that have one dipole on X pol. and another dipole (the
outlier probably) on Y pol.  Accepts three command line arguments:  
 1) stand number on X pol.
 2) stand number on Y pol.
 3) DRX filename(s)
At the end of the correlation, the visibilities are written to a NPZ file.

fringeBeam.py
-------------
Similar to fringeDipole.py, but expects the beam to be on X pol. and the dipole on
Y pol.  Accepts two command line arguments:
 1) stand number on Y pol.
 2) DRX filename(s)

fringeBeamHDF.py
----------------
Similar to fringeBeam.py, but the output is written to a single HDF5 file.  Accepts
two command line arguments:
 1) stand number on Y pol.
 2) DRX filename(s)

plotFringes.py
--------------
Simple script to load in a collection of NPZ files generated by fringeDipole.py/
fringeBeam.py and plot the visibility amplitude over time.

plotFringesHDF.py
-----------------
Similar to plotFringes.py, but reads in a HDF5 file created by fringeBeamHDF.py.

plotFringes2.py
---------------
Similar to plotFringes.py but displays a waterfall of the visibility and includes the
integrated bandpass.

plotFringes2HDF.py
------------------
Similar to plotFringes2.py, but reads in a HDF5 file created by fringeBeamHDF.py.

simulateFringesBright.py
------------------------
Given a collection of NPZ files generated by fringeDipole.py, simulate the fringes using 
the lsl.sim.vis.buildSimData() function and the bright sources listed in lsl.sim.vis.srcs.

readBinaryGainFile.py
---------------------
Simple script to parse a MCS/DP binary gain file (.gf) and print the gain matrix for each 
stand.

readBinaryDelayFile.py
----------------------
Simple script to parse a MCS/DP binary delay file (.df) and print the delays for each stand
in ns.

