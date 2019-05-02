.. _headreco_docs:

headreco
=========

Description
------------

headreco reconstructs a tetrahedral head mesh from T1- and T2-weighted structural MR images. It runs also with only a T1w image, but it will achieve more reliable skull segmentations when a T2w image is supplied.

.. attention:: headreco requires a MATLAB version compatible with SPM12 and CAT12. Linux and MacOSX users might need to configure MATLAB for usage with headreco, please see the :ref:`installation page <matlab_headreco>` page for more information. SPM12 and CAT12 are distributed together with SimNIBS.

Usage example
--------------

1. Open a terminal and go to the directory of the “Ernie” example data set.
2. Run the reconstruction:

.. code-block:: text

   headreco all --cat ernie org/ernie_T1.nii.gz org/ernie_T2.nii.gz

\
  The argument :code:`all` tells headreco to run all reconstruction steps including volume meshing. The subject ID (subID) :code:`ernie` is given next. Headreco will create a mesh named :file:`ernie.msh`, and a folder :file:`m2m_ernie/` that contains the segmentation results and the files needed for volume meshing. The input images are given as final arguments (first the T1, then the T2). The argument :code:`--cat` tells :code:`headreco` to use CAT12 for the segmentation of the WM and GM surfaces. The T2 image and CAT12 are optional, but highly recommended.

   Alternatively, the reconstruction can be run with only the T1w image as input, but this will result in a less accurate skull region:

.. code-block:: text

   headreco all --cat ernie org/ernie_T1.nii.gz

\

3. Check the results:

.. code-block:: text

   headreco check ernie

\
   This will show the reconstructed surfaces overlaid over the MR images using freeview. A second freeview will show the subject T1 registered to the  MNI template for visual inspection of the accuracy of the registration. In addition, you should have a look at the tetrahedral head mesh by loading it into Gmsh. In case freeview is not available, the spm viewer will be opened to allow for a basic check of the results.

Further notes
--------------

* When setting :code:`-d no-conform`, :code:`headreco` will not resample the supplied T1 and keep its original qform.
* Mesh resolution can be controlled using the -v option, which allows setting the vertex density (nodes per mm²) of the surface meshes. As a standard, :code:`headreco` uses 0.5 nodes per mm², resulting in head meshes with around 4 million tetrahedra.
* After the head mesh creation, temporary files are deleted to save diskspace. Adding :code:`--noclean` prevents this.
* Manual editing: Edit one or more of the binary masks stored in :file:`m2m_{subID}/mask_prep/`. Then run :code:`headreco surfacemesh subID` and :code:`headreco volumemesh subID` to re-create the head mesh based on the edited masks. Add :code:`--cat` to the surfacemesh step in case you used CAT12. Note: When using CAT12, surfaces instead of voxel masks will be stored for GM and WM in the :file:`mask_prep/` folder. For now, these surfaces cannot be manually improved.
* Transformation from and to MNI space: Both positions and results such as the electric field can be transformed between MNI and subject space. Please see below for a description of the corresponding command line programs. The transformation is based on a non-linear whole-head registration of the T1 of the subject to the MNI template that is determined during the SPM12 segmentation procedure. The transformations are stored in the :file:`m2m_{subID}/toMNI/` directory. Subject space is defined by the qform set in the :file:`m2m_{subID}/{subID}_T1fs_conform.nii.gz`, which can be found in the same folder as the head mesh. 
* When something goes wrong, you can check the :file:`m2m_{subID}/headreco_log.html` file.

References
-----------

`Nielsen, J. D., Madsen, K. H., Puonti, O., Siebner, H. R., Bauer, C., Madsen, C. G., ..., and Thielscher, A. (2018). Automatic skull segmentation from MR images for realistic volume conductor models of the head: Assessment of the state-of-the-art. NeuroImage, 174, 587-598. <https://doi.org/10.1016/j.neuroimage.2018.03.001>`_

