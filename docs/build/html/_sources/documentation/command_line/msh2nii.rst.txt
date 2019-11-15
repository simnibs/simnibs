.. _msh2nii_doc:

msh2nii
==========

Description
------------

Transforms fields in mesh files to nifti files.

Usage example
--------------

1. Open a simulation results folder, for example *ernie/simu* if you ran *ernie_simu.mat*

2. Run

.. code-block:: bash

  msh2nii ernie_TDCS_1_scalar.msh ../m2m_ernie/T1fs_conform.nii.gz ernie_TDCS_1

\

3. This will create the *ernie_TDCS_1_E.nii.gz*, *ernie_TDCS_1_J.nii.gz*, *ernie_TDCS_1_normE.nii.gz*, and *ernie_TDCS_1_normJ.nii.gz* files.
4. The files are in the same space as the *../m2m_ernie/T1fs_conform.nii.gz* file

Further notes
---------------

* Type :code:`msh2nii -h` for more information and options.
* This tool equivalent to selecting the **Interpolate to a nifti volume** option in the :ref:`GUI <sim_opt>` or set **map_to_vol** to *true* in the :ref:`SESSION structure <session_doc>`.
* Can also be used to create tissue masks with the :code:`--create_masks` argument.
