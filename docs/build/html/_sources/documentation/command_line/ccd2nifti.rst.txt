.. _ccd2nifti_doc:

ccd2nifti
==========

Description
------------

Kris will tell you what it is

Usage example
-------------

!! This needs update from Kris!!
1. Rigid registration: open a terminaĺ and run

.. code-block:: bash

  register -f ernie_T1.nii.gz -m ernie_T2.nii.gz -o ernie_T2_registered.nii.gz -dof 6


2. Affine registration: open a terminal and run

.. code-block:: bash

   register -f ernie_T1.nii.gz -m MNI152_T1_1mm.nii.gz -o MNI_template_registered.nii.gz -dof 12


Further notes
---------------

* Type :code:`ccd2nifti -h` for the command line help
