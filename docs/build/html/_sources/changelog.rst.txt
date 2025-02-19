.. _changelog:

Changelog
===========
4.5.0 
------
  * New optimization method TMS that also supports bent and flexible coils, thereby systematically avoiding intersections of the coil with the head. 
  * New format (.tcd) for TMS coils that supports flexible and multi-element coils. Coils can now also be defined by their winding geometries. The A-fields (magnetic vector potential) are then determined via numerical integration of line integrals. This method complements coil definitions using magnetic dipoles and using precalculated A-fields stored on a regular grid (i.e. as NIfTI). The tcd-format supports all three cases.
  * TMS coil models for Brainsway H1, H4 and H7, and for MagVenture MST twin coil.
  * New optimization method for TES, including montages with rectangular electrodes, center-surround montages, temporal interference stimulation and electrode arrays for tumor treating field therapies.
  * New class to support the convenient definition of regions-of-interest.
  * New dataset Ernie Extended.
  * New non-human primate dataset.
  * Tutorial for calculating :ref:`EEG leadfields <eeg_leadfields>` with SimNIBS for use in `FieldTrip <https://www.fieldtriptoolbox.org/>`_ and `MNE-Python <https://mne.tools/stable/index.html>`_ added.
  * Added `MUMPS solver <https://mumps-solver.org/index.php>`_ as option for the FEM calculations on Macs with Apple Silicon
  * Added custom compiled version of `FMM3d <https://github.com/flatironinstitute/FMM3D>`_ so that it runs on Macs with Apple Silicon
  * Included JupyterLab (start via ``simnibs_jupyter`` on command line)
  * Update to python 3.11 and corresponding updates of most of the included packages.

NOTES & Known issues:
 * The PARDISO solver does not work on Apple Silicon (use MUMPS instead).
 * Support of SimNIBS on MacOS on Intel is discontinued.
 * Installation fails on paths with non-standard characters, such as backslash, chinese characters, ... (workaround: provide another path)
 * simnibs_gui does not start on some linux systems, e.g. with wayland (workaround: ld preloading of libstdc++.so.6 seems to help; example: export LD_PRELOAD=/usr/lib/libstdc++.so.6 - path needs to be adjusted according to library path on local system)
 * Setting the number of cpus in run_simnibs will be only taken into account in the new tms_flex_opt and tes_flex_opt optimizations, but not the previous TMS optimization methods.
 * Setting custom conductivities and defining custom tissue types in the new tms_flex_opt and tes_flex_opt does not work yet.


4.1.0 
------
 * Tetrahedral quality of the meshes was increased substantially to improve numerical accuracy of the FEM calculations and remove outliers in the calculated electric fields
 * Option *--fs-dir* added to charm to use white matter and pial surfaces from FreeSurfer for more accurate representation of smaller sulci in the head meshes
 * I/O functions for neuronavigation data have been updated to support new Brainsight version 2.5.3
 
4.0.1
------
  * changed Brainsight position import/export to support only NIfTI:Aligned (to avoid ambiguities)
  * small bug fixes

4.0.0
------
 * New head segmentation and meshing pipeline *charm* with improved accuracy and robustness
 * Many new TMS coil models
 * New head models with additional tissue types, in particular spongy bone and large blood vessels
 * New flexible meshing approach to simplify manual editing and inclusion of custom tissue types in the head mesh
 * New command line tool *meshmesh* to support meshing of custom geometries
 * Support of Nx1 center-surround montages in python and Matlab
 * Added python examples for temporal interference simulations using pre-calculated leadfields for speed
 * Added basic I/O functions for neuronavigation data in python
 * Update to python 3.9
 * Major code cleanup and restructuring under the hood
 * Tested on Windows 10, Linux and Macs with Intel and Apple Silicon
 * Headreco and mri2mesh are deprecated.
 
NOTE: Simnibs 4 is NOT backwards compatible. Head models created with charm cannot be used in older versions. Head models from older versions need to be converted using the command line tool *convert_3_to_4* for use with SimNIBS 4.
 
Known issues:
 * Installation fails on paths with non-standard characters, such as backslash, chinese characters, ... (workaround: provide another path)
 * simnibs_gui does not start on some linux systems, e.g. with wayland (workaround: ld preloading of libstdc++.so.6 seems to help; example: export LD_PRELOAD=/usr/lib/libstdc++.so.6 - path needs to be adjusted according to library path on local system)
 * The PARDISO solver (available as option for the FEM calculations) does not work on Apple Silicon
 
3.2.6
------
* fixed a bug causing headreco to fail on MacOS
* Note: On Apple M1s using .ccd coil files does not yet work, please use niftis instead.

3.2.5
------
* added skin smoothing options for TMS optimization
* more informative matlab error messages
* small bug fixes

3.2.4
------
 * Small bug fix related to gmsh options.

3.2.3
------
 * Gmsh version changed to avoid issues with Big Sur (only MacOS)
 * Headreco bug fixes to make meshing more stable in the eye region and air cavities
 * Added a Nx1 example for a center-surround electrode montage

3.2.2
------
 * Added matlab and python examples for TACSchallenge
 * Added matlab examples to calculate TI amplitudes
 * TMS optimization now prints optimal position in log file
 * Added reading of leadfields in mesh_load_hdr5.m
 * Added work-around to enable installation on Mac OS 11 Big Sur; NOTE: If you attempted to install a previous version of SimNIBS 3.2, you have to wipe the installation folder (in "/Users/username/Applications/SimNIBS-3.2" before attempting to install again; NOTE: This work-around is temporary and will be removed once the issues in the underlying python packages have been resolved
 * Changed the skin smoothing iterations in headreco from 5 to 20, which should result in a smoother skin surface. Note: this changes the standard behavior of headreco, so results might differ slightly.
 * Per default, headreco does not print the cat summary pdf anymore. The --cat_no_print flag was removed; Instead, use --cat_print now in case you need the summary.


Bug fixes:
 * TDCS Network Optimization: Fixed bugs to accept images with NaNs, binary images, and images of size NxMxKx1; weights of eyes are now set to 0.
 * Added buffered read for gmsh v2 files in python to resolve speed issue on clusters
 * Resolved a bug causing some points of the individual middle gm surface to be falsely interpreted as outside gm when interpolating results to individual gm surface
 * Further small fixes across the code
 * Fixed electrode meshing that sometimes caused some parts of the electrodes to be detached.

 
Known issues:
 * mri2mesh does not work with Freesurfer 7; please use Freesurfer 6 for now
 * SimNIBS is so far not tested on Macs with Apple Silicon, and is likely to give errors on those machines


3.2
----
 * Added Auxiliary Dipole Method (ADM) TMS optimization (contributed by Luis Gomez)
 * Added TES magnetic field calculations for MRCDI/MREIT (contributed by Hassan Yazdanian)
 * Added TES optimization with field strength
 * Added TES optimization for brain network targeting
 * FMM-based coil A field calculations from :file:`.ccd` files
 * Refactoring of optimization code
 * Calculating coil-cortex distances during TMS simulations
 * New installers


3.1
----
 * Added TDCS and TMS Optimization
 * Added support to the MKL PARDISO solver
 * Minor bug fixes


3.0
-----
 * Major update to SimNIBS
 * New integrated solver based on PETSc and Hypre with huge speed ups!
 * New installation procedure
 * Changed headreco call (:code:`-d no-conform` and CAT12 now standard)
 * New coil models
 * Moved to Python 3.7
 * Updated documentation
 * Uncertainty quantification support
 * Improved results visualizations
 * SimNIBS is now installable as a python package


2.1.2 
---------
  * This upgrade focused on the MATLAB library for SimNIBS
  * The MATLAB library can be found in in the *matlab/functions/* folder in the simnibs directory
  * The MATLAB examples have been expanded to show the new features
  * We also made a few bug fixes to *get_fields_at_coordinates* and the graphical user interface

2.1.1
---------
  * This upgrade focused on usability of SimNIBS
  * Simnibs is now shipped with spm12 and cat12. Installing those separately is no longer necessary to run **heareco**. When you have Matlab installed, you're ready to go after the SimNIBS installation.
  * Added features for automatically calculate eeg positions and easily using them to set-up simulations in the GUI. Added a script called **get_eeg_positions** in order to use the new features on head models ran with SimNIBS 2.1.0
  * Improved Matlab and python scripting. Examples can be fount in $SIMNIBSDIR/matlab/examples and $SIMNIBSDIR/python_examples
  * Changed the sign of the normals when interpolating to the cortical surface
  * the SimNIBS python installation no longer requires scikit-image

2.1.0
---------
  * 2.1.0 is a major update of SimNIBS 2
  * New head segmentation script **headreco**
  * New post-processing options to transform fields to NIfTI volumes, MNI space, FreeSurfer overlays and FsAverage space
  * New scripts to calculate EEG 10-10 positions
  * New MATLAB library, including example scripts e.g. to set up simulations for ring electrodes
  * New example data sets, including an extended MNI template
  * Major refactoring under the hood, for a cleaner experience and quicker future updates
  * **Head segmentations and simulation files created with SimNIBS 2.0 are incompatible with SimNIBS 2.1**

2.0.1g
---------
  * More fixes to the GUI
  * Changes address to Miniconda during installation procedure
  * Changed bug in simnibs.py where it would look for files that didn't exist

2.0.1f
---------
  * Changed 3dcalc wrapper in Linux
  * Fixed bugs in the GUI related to the PySide->PyQt changes

2.0.1e
---------
  * Changed from PySide to PyQt
  * Changed getopt in osx to a wrapper script which will call getopt_o adjusting DYLD_LIBRARY_PATH
  * Changed scalp color on GUI

2.0.1d
---------
  * Fixed the intallation of qt on mac
  * Fixed an installation bug that occurred when reinstalling simnibs without starting a new terminal window

2.0.1c
---------
  * freeglut is no longer required
  * now a local verision of libXp, libXpm and libXmu is provided. This should make the installation easier
  * solved bug in the anisotropic conductivity calculations that would happen if there was a "." in the path


2.0.1b
---------
  * Now SimNIBS gui gives out a warning if there are any spaces in the file path
  * The installation procedure will now install freeglut on Linux

2.0.1a
---------
  * support of conductivity tensors for gray and white matter added to GUI
  * script dwi2cond added to estimate conductiviy tensors from diffusion MRI
  * automatic installation procedure changed to use miniconda
  * bug in mri2mesh fixed which prevented it to use the T2 image to reconstruct the skull
  * Changed standard colors in GUI
  * The GUI now lets you set TMS coil distances
  * Removed deprecated post processing options 
  * Added new TMS coil files
  * Fixed bug where the electrode thickness would change every time the electrode edition window opened
  * Fixed bug where every simulation would use the same conductivities
  * GUI now supports advanced electrode modeling
  * **ATTENTION** old .simnibs files maybe incompatible with the new version

