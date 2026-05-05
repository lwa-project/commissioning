Delay Calibration Guide
=======================

Before You Start
----------------
 * Make sure you have an SSMIF with a priori cable delays measured via TDR
 * Have an initial mask of bad antennas in the SSMIF
 * Make sure you can capture *good* TBT data from *all* antennas

Pass 0
------
 * If you have TBT data dumped locally, combine the individual files into a
   single file with `tbxMux.py`
 * Choose a calibration strategy of either a point source model-based calibration
   using Cyg A or a diffuse emission-based calibration using observations when
   the Galactic plane is low.  `classifyData.py` can be used to help determine
   which data are useful for which approach.
 * Establish a baseline SEFD value using `estimateSEFD.py` and `fitEstimatedSEFD.py`
   with data near transit of Cyg A

Pass 1
------
 * Correlate the data using LSL's `correlateTBX.py` function with your starting
   SSMIF file and the `-2` option so that you get both XX and YY.
 * Run either `applySelfCalTBX.py` (Cyg A) or `applyDiffuseSelfCalTBX.py` (diffuse
   emission) on the data.  These scripts will generate `.sc` delay calibration
   files that contain per-stand and per-polarization delays in nanoseconds.
 * Combine all of the `.sc` files into a single `merged.delays` file using
   `checkConsistency.py`.
 * Use `estimateSEFD.py` and `fitEstimatedSEFD.py` with the `-c` option to
   evaluate how well the merged delays apply to the data.
    * If the SEFD and pointing error improve continue on.
    * Otherwise, try running the calibration script with additional options
      (smaller frequency range, auto-correlation-based filtering with `-f`,
      regularization with `-e`).
    * Also try building the merged delays with only `.sc` files that have
      fully converged (check the `Converged XX`/`Converged YY` lines in the
      `.sc` file headers).
 * Convert the delays into stretch factors suitable for updating the starting
   SSMIF using `convertDelayToStretch.py`
 * Use the starting SSMIF and stretch factors to generate a new SSMIF with
   `applyNewStretchFactors.py`.  Also, consider using any stand flags generated
   by the calibration script to update the bad antenna mask in the new SSMIF.

Passes 2...N
------------
 * Correlate the data again using the updated SSMIF from the previous pass.
   From pass 3 onward, add `-d 4` to `correlateTBX.py` to decimate the
   correlator output by a factor of 4 in frequency.  The wider self-cal
   bands on later passes don't need full channel resolution to fit delays,
   and decimation cuts FITS IDI size and downstream processing time.
 * Run the calibration on the data again.  If you used a smaller frequency
   range on the previous pass consider trying a larger range now.  If you
   used regularization consider trying less regularization now.  etc.
 * Combine the `.sc` files into a `merged.delays` file.
 * Evaluate the SEFD and pointing error.
 * Convert to stretch factors and generate another updated SSMIF.
 * Continue iterating until you've reached a reasonable sensitivity.  ~12 kJy
   for a full station and ~50 kJy for a 64-element mini-station.

Example `applySelfCalTBX.py` Flags
-----------------------------------
A rough flow of iterations through `applySelfCalTBX.py` might be:
 1. First pass with 65-70 MHz, filtering and regularization: `-l 65 -u 70 -f -e 0.5`
 2. Second with a wider frequency range: `-l 65 -u 75 -f -e 0.5`
 3. Third even wider: `-l 55 -u 75 -f -e 0.5`
 4. Fourth even wider: `-l 35 -u 75 -f -e 0.5`
 5. Fifth with less regularization: `-l 35 -u 75 -f -e 1`
