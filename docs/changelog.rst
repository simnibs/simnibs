Changelog
===========

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

