.. _nnav_docs:

Neuronavigation Data Import/Export
==================================

.. toctree::
   :maxdepth: 1
   
   localite
   brainsight
   softaxic

General Notes
-------------
- The import from SimNIBS coil positions/orientations to neuronavigation softwares and vice versa only works if the same T1 has been used.
- In addition, SimNIBS expects the :code:`qform` and :code:`sform` transformation information in the NIfTI header to be equal. If you used the :code:`--forceqform` flag during the meshing, make sure to use the SimNIBS processed :code:`T1.nii.gz` image for the neuronavigation. Otherwise a correct transformation of coil positions/orientations cannot be guaranteed!