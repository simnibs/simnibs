.. _eeg_positions_doc:

eeg_positions
===============

Description
-------------

Calculates the 10/10 EEG positions based on 4 fiducial points based on the UI 10/10 definition provided in (`Jurcak et al., 2007 <https://doi.org/10.1016/j.neuroimage.2006.09.024>`_).⁠

Usage example
--------------

1. Open a terminal and go to the directory of the “Ernie” example data set.
2. Open the mesh file in the GUI

.. code-block:: bash

   simnibs_gui ernie.msh

\ 

3. Visually determine the positions of the nazion (Nz), left preauricular point (LPA), right preauricular point (RPA) and inion (Iz) on the skin surface and click on them to get their positions. Write down the coordinates displayed in the bottom of the window.

4. Run using the coordinates taken in step 3

.. code-block:: bash

  eeg_positions -m ernie.msh -o ernie_10_10 -Nz 1.6 97.1 3.4 -LPA -82.2 -14.3 -22.2 -Iz -4.8 -108.8 -32.7 -RPA 79.7 -18.9 -24.6

\

5. Two files will be created: :file:`ernie_10_10.csv` and :file:`ernie_10_10.geo`

6. Check the positions by opening the mesh file in Gmsh, going to *File* →  *Merge* and opening the :file:`ernie_10_10.geo` file for a visual check of the electrode positions

7. :file:`ernie_10_10.csv` is a text file containing names and positions of the EEG electrodes. For now, you have to manually enter these coordinates when defining electrodes in simnibs_gui. Alternatively, you can write your own matlab- or python-script that uses the positions to create a .mat-file with the simulation information.

Further notes
--------------

* We have compared positions obtained using this procedure against positions automatically obtained from MNI electrodes definitions and found no large differences (`Saturnino et al., 2018 <https://doi.org/10.1101/500314>`_). The later are automatically calculated during our head meshing programs **headreco** and **mri2mesh** and stored in the :file:`m2m_{subID}/eeg_positions/EEG10-10_UI_Jurak_2007.csv` file.


