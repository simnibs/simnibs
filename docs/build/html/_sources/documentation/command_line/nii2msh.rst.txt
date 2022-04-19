.. _nii2msh:

nii2msh
==========

Description
------------

Interpolates NifTI volumes to mesh files. 

Usage example
-------------

1. Open a terminal and go to the directory of the “Ernie” example data set.
2. Run

.. code-block:: bash

  nii2msh m2m_ernie/T1.nii.gz m2m_ernie/ernie.msh ernie_with_T1.msh

\

3. The *ernie_with_T1.msh file* will contain the intensity values of the T1 image as data field.

Further notes
---------------

* Type :code:`nii2msh -h` for more information and options
* The -ev option will map the eigenvectors associated with the largest eigenvalues of a tensor to the mesh, which can be useful for visual control of the vector orientations (e.g. :code:`nii2msh -ev m2m_ernie/DTI_coregT1_tensor.nii.gz m2m_ernie/ernie.msh ernie_tensor.msh` adds the first conductivity eigenvector to the mesh).
* The data will be added to the elements (tetrahedra) as default. To add it to the nodes, use the option :code:`--type nodes`.



