.. _mni2subject_coords_docs:

mni2subject_coords
===================

Description
------------

Transforms coordinates from MNI space to subject space.

Usage example
---------------


.. code-block:: bash

 mni2subject_coords -m m2m_ernie/ -s list.csv -o ernie_10_10.csv

\
  This will transform a list of defined in MNI space saved in the :file:`list.csv` file to subject space.

Further notes
---------------

* Type :code:`mni2subject_coords -h` for more information and options, including on the format of the *.csv* files.
* There are many different behaviors of this script depending on specifications of the *.csv* file. Type :code:`mni2subject_coords -h` for a full description
* This script is used to determine the electrode positions in the *m2m_{subID}/eeg_positions/EEG10-10_UI_Jurak_2007.csv*
* When transforming electrode positions, the results are always projected to the skin after the transformation.


