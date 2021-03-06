DRX Commissioning Scripts
=========================

decodeStatus.py
---------------
Take the status code returned by a DP_ BOARD_STAT query and break it down into its
constituent parts.

gain.py
-------
JPL-created script for creating gain files to be used by DRX.

delay.py
--------
JPL-create script for created delay files to be used by DRX.

test_analog_in2.py
------------------
Script for setting up a test DRX observation with a variety of tuning, frequency, 
bandwidth, and gain settings.  Needs gain.py and delay.py to work.

generateDelays.py
-----------------
Generate gain and delay files for a particular topocentric pointing (azimuth and 
elevation in degrees) at a specified frequency.

.. note::
	genernateDelays.py generates text files for gains and delays that need to
	be converted to binary packed files with the MCS/Executive commands megfg 
	and medfg.

estimateBeam.py
---------------
Use the LSL misc.beamformer's integer delay and sum to estimate the beam shape 
towards a particular azimuth and elevation and save the beam to a Numpy NPZ file.

driftcurve.py
-------------
Using the beam created by estimateBeam.py, generate a temperature as a function of 
time for that pointing over an entire day.  This script is based on the LSL 
driftcurve.py script.

checkTimetags.py
----------------
Simple script to check the time tag sequence in a DRX file and look for strange 
jumps and other badness.  It can also be used to check the DRX decimate (sample 
rate) field and if the time has been set correctly on the various boards.

fastDRXCheck.py
---------------
Script to quickly check a DRX file for time tag problems and, optionally, split the
file at any time tag problems.

getTuningOffset.py
-------------------
Check for an offset in the time tags associated with the two tunings beyond what is
controlled by the time offset stored in the frame header.  This should read in the 
first complete set of frame from all four beamtunepols and compare tunings 1 and 2.  
Any differences between the times are noted.

checkTuningCoherency.py
-----------------------
Given a a DRX file with both tunings set to the same parameters, check for coherency
across the tunings via cross-correlation.  For each set of data, print out the lag
(in DRX samples) where the cross-correlation function has its peak and the normalized
correlation magnitude relative to the previous set of tuning 1 data.

beamCoherency.py
----------------
Given two or more DRX files from different beams, check for coherency between the 
beams and make sure that the beams agree with the T_NOM values.  For each 
beam/tuning/pair, print out the lag (converted to ticks), the time tag differences
(both raw and T_NOM corrected) and the normalized correlation magnitude relative to 
the same-beam tuning cross-correlation.

drxTimeseries.py
----------------
Script for plotting time series data from a DRX file.  This file is also available 
as part of the LSL package (version >= 0.4.0) but is included here since this is 
the script we actually used at JPL for pre-acceptance testing.

drxPower.py
-----------
Modified version of drxTimeseries.py that plots the power for each DRX sample instead
of the raw I/Q values.  The power is summed for the integration time (listed as 
average) and displayed for the duration that it is computed.

drxFileCheck.py
---------------
Script for inspecting one second of DRX data every 15 minutes to check the gain
settings.

eyeDiagram.py
-------------
Create an eye diagram for some portion of a DRX file.

eyeDiagram2.py
--------------
Look for glitches in a DRX file by fitting a quantized sine wave to the data.

Fringing
--------
See the Fringing/README file.

HDF5
----
See the HDF5/README file.

