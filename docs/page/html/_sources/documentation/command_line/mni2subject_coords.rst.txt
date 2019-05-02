.. _mni2subject_coords_docs:

mni2subject_coords
===================

Description
------------

Transforms coordinates from MNI space to subject space.

Usage example
---------------

1. Open a terminal and go to the directory of the “Ernie” example data set.
2. Run

.. code-block:: bash

 mni2subject_coords -m m2m_ernie/ -s $SIMNIBSDIR/resources/ElectrodeCaps_MNI/EEG10-10_UI_Jurak_2007.csv -o ernie_10_10.csv

\
  This will transform a list of electrode positions that were defined in MNI space according to the EEG 10/10 system (using the “UI” convention as defined in `Jurcak et al., 2007 <https://doi.org/10.1016/j.neuroimage.2006.09.024>`_) to subject space. The transformed positions are stored in the file *ernie_10_10.csv*.

Further notes
---------------

* Type :code:`mni2subject_coords -h` for more information and options, including on the format of the *.csv* files.
* There are many different behaviors of this script depending on specifications of the *.csv* file
* This script is used to determine the electrode positions in the *m2m_{subID}/eeg_positions/EEG10-10_UI_Jurak_2007.csv*
* When transforming electrode positions, the results are always projected to the skin after the transformation.


