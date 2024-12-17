.. _charm_docs:

charm
=========

.. note:: Charm replaces :ref:`headreco_docs` or :ref:`mri2mesh_docs` that were included in older SimNIBS versions.

Description
------------

charm reconstructs a tetrahedral head mesh from T1- and T2-weighted structural MR images. It runs also with only a T1w image, but it will achieve more reliable skull segmentations when a T2w image is supplied.

Usage example
--------------

1. Open a terminal and go to the :file:`ernie/` folder of the example data set.
2. Run the reconstruction:

  .. code-block:: bash

     charm ernie org/ernie_T1.nii.gz org/ernie_T2.nii.gz

  \
  The subject ID (subID) :code:`ernie` is given as first argument. Charm will create a folder named :file:`m2m_ernie/` that contains the segmentation results and the final head mesh :file:`ernie.msh`. The input images are given as final arguments (first the T1, then the T2).

\

  Alternatively, the reconstruction can be run with only the T1w image as input, but this can result in a less accurate skull region:

  .. code-block:: bash

     charm ernie org/ernie_T1.nii.gz

  \

3. Check the segmentation. Click on the final segmentation viewer in the results.html (to be found in the m2m-folder of the subject). The viewer shows the outlines of the reconstructed tissue compartments, enabling a visual check whether the outlines are accurate.

Further notes
--------------

* If you encounter spurious segmentation results this *could* be due to a suboptimal affine registration between the anatomical image(s) and the atlas. Please see the tutorial :ref:`fix_affine_registration_tutorial`.
* charm can use the cortical surfaces created by FreeSurfer recon-all to achieve a more accurate representation of smaller sulci in the head meshes (option *--fs-dir*).						
* Please see the tutorial :ref:`fixheadmodel_tutorial` in case manually fixes to the segmentation are needed.


References
-----------

`Puonti O, Van Leemput K, Saturnino GB, Siebner HR, Madsen KH, Thielscher A. (2020). Accurate and robust whole-head segmentation from magnetic resonance images for individualized head modeling. Neuroimage, 219:117044. <https://doi.org/10.1016/j.neuroimage.2020.117044>`_
