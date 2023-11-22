.. _coords_doc:

Coordinates in SimNIBS
========================


Definitions
------------

Coordinates in SimNIBS are **always world coordinates in subject space**. Coordinates are given in **millimiters** and along *x*, *y* and *z* axis. The origin and the alignment of the axes correspond to those of the original T1w image that was used for head model creation in *charm* (stored as :file:`m2m_{subID}/T1.nii.gz`).


MNI Transformations
----------------------

During head meshing, SimNIBS calculates 6 degrees-of-freedom, 12 degrees-of-freedom and non-linear MNI transformations. They are stored in the :file:`m2m_{subID}/toMNI/` folder.
SimNIBS uses these transformations when transforming simulation results to MNI space (see the :ref:`sim_opt` and :ref:`map_to_mni attribute <session_doc>`).

Command line Utilities
~~~~~~~~~~~~~~~~~~~~~~~

SimNIBS has several command line utilities to transform data and coordinates between MNI and Subject space

* :ref:`mni2subject<mni2subject_docs>`, to transform data in MNI space to subject space.
* :ref:`subject2mni <subject2mni_docs>`, to transform data in subject space to MNI space.
* :ref:`mni2subject_coords <mni2subject_coords_docs>`, to transform coordinates from MNI space to subject space.
* :ref:`subject2mni_coords <subject2mni_coords_docs>`, to transform coordinates from subject space to MNI space.

Python and MATLAB Functionas
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
We also provide *Python* and *MATLAB* interfaces to transform coordinates between MNI and Subject Space. Those have the signature

* :code:`mni2subject_coords(coords_mni, subdir, transformation_type)`

  * :code:`coords_mni` is a set of coordinates in MNI space
  * :code:`subdir` is the subject segmentation directory, for example :file:`m2m_ernie`
  * (Optional) :code:`transformation_type` is the type of transformation, can be :code:`nonl` (default), :code:`12dof` or :code:`6dof`.

And equivalently to :code:`subject2mni_coords`.
