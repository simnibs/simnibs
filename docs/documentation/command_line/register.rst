.. _register_doc:

register
==========

Description
------------

A lightweight tool to perform rigid (--dof 6) or affine (--dof 12) registration. Outputs a resampled image and a the world-to-world transformation matrix in RAS space. *NOTE* The rigid registration supports inter-modality (across contrast) registration while the affine currently only supports intra-modality registration (the same contrast).

Usage example
-------------
1. Rigid registration: open a terminaĺ and run

.. code-block:: bash

  register -f ernie_T1.nii.gz -m ernie_T2.nii.gz -o ernie_T2_registered.nii.gz -dof 6


2. Affine registration: open a terminal and run

.. code-block:: bash

   register -f ernie_T1.nii.gz -m MNI152_T1_1mm.nii.gz -o MNI_template_registered.nii.gz -dof 12


Further notes
---------------

* Type :code:`register -h` for more information and options
