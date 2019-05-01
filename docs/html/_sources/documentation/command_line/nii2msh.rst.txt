.. _nii2msh:

nii2msh
==========

Description
------------

Interpolates NifTI volumes to mesh files. Mostly used internally by **dwi2cond**

Usage example
-------------

1. Open a terminal and go to the directory of the “Ernie” example data set.
2. Run

.. code-block:: bash

  nii2msh -ev d2c_ernie/dti_results_T1space/DTI_conf_tensor.nii.gz ernie.msh ernie_tensor.msh

\

3. This will create the *ernie_tensor.msh file*. This file will contain the eigenvectors associated with the largest eigenvalues of the conductivity tensors (for visual control).

Further notes
---------------

* Type :code:`nii2msh -h` for more information and options

