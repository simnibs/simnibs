.. _coords_doc:

Coordinates in SimNIBS
========================


Definitions
------------

Coordinates in SimNIBS are **always world coordinates in subject space**. Coordinates are given in **millimiters** and along *x*, *y* and *z* axis.

The origin and the alignment of the axes can vary in head models ran with :ref:`headreco <headreco_docs>`. In head models ran with :ref:`mri2mesh <mri2mesh_docs>`, the coordinates are based on the conform coordinate system.


Voxel to SimNIBS Coordinates
-----------------------------

The head meshing procedures create a NifTI file called :file:`m2m_{subID}/T1fs_conform.nii.gz`. This file contains the input T1 image of the subject with a changed header such that the world coordinates reflect those used by SimNIBS.

In :ref:`headreco <headreco_docs>`, you can avoid re-sampling to conformed space by setting the :code:`-d no-conform` option. This way, SimNIBS will use the same coordinate system as the input NifTI file.


MNI Transformations
----------------------

During head meshing, SimNIBS calculates, 6 degrees-of-freedom, 12 degrees-of-freedom and non-linear MNI transformations. They are stored in the :file:`m2m_{subID}/toMNI/` folder.
SimNIBS uses these transformations when transforming simulation results to MNI space (see the :ref:`sim_opt` and :ref:`map_to_mni attribute <session_doc>`) as well as in the command line utilities

* :ref:`mni2subject<mni2subject_docs>`, to transform data in MNI space to subject space.
* :ref:`subject2mni <subject2mni_docs>`, to transform data in subject space to MNI space.
* :ref:`mni2subject_coords <mni2subject_coords_docs>`, to transform coordinates from MNI space to subject space.
* :ref:`subject2mni_coords <subject2mni_coords_docs>`, to transform coordinates from subject space to MNI space.

