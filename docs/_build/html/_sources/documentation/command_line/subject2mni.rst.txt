.. _subject2mni_docs:

subject2mni
============

Description
------------

Transforms the simulation results stored in a mesh or nifti file to MNI space.

Usage example
---------------

1. Open a simulation results folder, for example :file:`ernie/simu` if you ran :file:`ernie_simu.mat`

2. Run

.. code-block:: bash

  subject2mni -i ernie_TDCS_1_scalar.msh -m ../m2m_ernie/ -o ernie_TDCS_1

\

  This will create the :file:`ernie_TDCS_1_MNI_E.nii.gz`, :file:`ernie_TDCS_1_MNI_J.nii.gz`, :file:`ernie_TDCS_1_MNI_magnE.nii.gz`, and :file:`ernie_TDCS_1_MNI_magnJ.nii.gz` files, with the fields transformed to MNI space.

Further notes
---------------

* Type :code:`subject2mni -h` for more information and options 
* This tool equivalent to selecting the **Transform to MNI space** option in the :ref:`GUI <sim_opt>` or **map_to_mni** to *true* in the :ref:`SESSION <session_doc>` structure.


