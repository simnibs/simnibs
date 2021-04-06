.. _msh2nii_doc:

msh2nii
==========

Description
------------

Transforms fields in mesh files to nifti files.

Usage example
--------------

1. Open a simulation results folder. For the following, it is assumed that you ran a TDCS simulation on ernie.

2. Run

.. code-block:: bash

  msh2nii ernie_TDCS_1_scalar.msh ../m2m_ernie/T1.nii.gz ernie_TDCS_1

\

3. This will create one file for each data field in the mesh, e.g. *ernie_TDCS_1_E.nii.gz*, *ernie_TDCS_1_J.nii.gz*, *ernie_TDCS_1_magnE.nii.gz*, and *ernie_TDCS_1_magnJ.nii.gz* ...
4. The files are in the same coordinate space as the *../m2m_ernie/T1.nii.gz* file

Further notes
---------------

* Type :code:`msh2nii -h` for more information and options.
* This tool performs an equivalent function to selecting the **Interpolate to a nifti volume** option in the :ref:`GUI <sim_opt>` or set **map_to_vol** to *true* in the :ref:`SESSION structure <session_doc>`.
* Can also be used to create tissue masks with the :code:`--create_masks` argument (one binary mask per tissue).
* Using the :code:`--create_label` argument will create an image with the tissue labels in the mesh.
