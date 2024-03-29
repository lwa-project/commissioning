DRX Processing Scripts with HDF5 Output
=======================================

hdfWaterfall.py
---------------
Version of drxWaterfall.py that writes to HDF5 rather than NPZ to provide
an easier path to getting the data into other packages, e.g., Matlab.  This
script supports both linear polarization products and Stokes parameter.

drspecCheckTimetags.py
----------------------
Script to check the flow of time in a DR spectrometer data file.

drspecFileCheck.py
------------------
Script to scan a DR spectrometer file and look at the data quality.

drspec2hdf.py
-------------
Simple script to convert a binary DR spectrometer file into an HDF5 file that
is compatible with the above HDF5 organization.

splitObservations.py
--------------------
Read in a DRX/HDF5 watefall file and split out various observations.  The 
observations can be split by:
  * Target Source
  * Observation ID number
or the script can be used to just list the observations within an HDF5 file.

calculateSK.py
--------------
Read in a DRX/HDF5 wataerfall file and calculate the pseudo-spectral kurtosis 
(pSK) for each observation (XX and YY or I).  The pSK values are "pseudo" since 
the spectra contained within the HDF5 file are averaged over N FFTs where 
N > 1.  The pSK value for each channel is calculated by examining M consecutive 
spectra.

Note:  Although the pSK method works reasonable well for short integration times
(<~0.3 s) the method may break down for much longer integration times.

dedisperseHDF.py
----------------
Given an HDF5 file, apply incoherent dedispersion to the data at the specified 
dispersion measure and save the results to a new file.

decimateHDF.py
--------------
Given an HDF5 file, decimate the data contained in it in both time and frequency, 
and save the results to a new file.

linear2stokes.py
----------------
Given an HDF5 file containing linear polarization products, compute the related 
Stokes parameters and save the results to a new file.  The combinations supported
are:
  * XX and YY -> I and Q
  * XY_real and XY_imag -> U and V

stokes2linear.py
----------------
Given an HDF5 file containing Stokes parameters, compute the related linear 
polarization products and save the results to a new file.  The combinations supported
are:
  * I and Q -> XX and YY
  * U and V -> XY_real and XY_imag

plotHDF.py
----------
Version of plotWaterfall.py for viewing HDF5 files created by hdfWaterfall.py and
drspec2hdf.py.  This viewer supports all linear and Stokes data products.
