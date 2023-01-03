.. _headreco_docs:

headreco
=========

.. warning:: *Headreco* is deprecated. Please use :ref:`charm_docs` instead. 

Description
------------

headreco reconstructs a tetrahedral head mesh from T1- and T2-weighted structural MR images. It runs also with only a T1w image, but it will achieve more reliable skull segmentations when a T2w image is supplied.

.. attention:: headreco requires a MATLAB. Linux and MacOSX users might need to configure MATLAB for usage with headreco, please see the :ref:`see here <matlab_setup>` for more information.

.. note:: headreco relies on `SPM12 <https://www.fil.ion.ucl.ac.uk/spm/software/spm12/>`_ and `CAT12 <http://www.neuro.uni-jena.de/cat/>`_ for segmentation. Both libraries are distributed together with SimNIBS, so you don't need to download anything extra

Usage example
--------------

1. Open a terminal and go to the :file:`ernie/` folder of the example data set.
2. Run the reconstruction:

  .. code-block:: bash
  
     headreco all ernie org/ernie_T1.nii.gz org/ernie_T2.nii.gz
  
  \
  The argument :code:`all` tells headreco to run all reconstruction steps including volume meshing. The subject ID (subID) :code:`ernie` is given next. Headreco will create a mesh named :file:`ernie.msh`, and a folder :file:`m2m_ernie/` that contains the segmentation results and the files needed for volume meshing. The input images are given as final arguments (first the T1, then the T2).

\

  Alternatively, the reconstruction can be run with only the T1w image as input, but this will result in a less accurate skull region:

  .. code-block:: bash
  
     headreco all ernie org/ernie_T1.nii.gz
  
  \

  You can also run headreco **without** using CAT12 by calling.

  .. code-block:: bash
  
     headreco all --no-cat ernie org/ernie_T1.nii.gz
  
  \

  Segmentations without CAT12 will run faster, but will have a poorer reconstruction of the cortical gray and white matter.

3. Check the results:

  .. code-block:: text
  
     headreco check ernie
  
  \
   This will show the reconstructed surfaces overlaid over the MR images using freeview. A second freeview will show the subject T1 registered to the  MNI template for visual inspection of the accuracy of the registration. In addition, you should have a look at the tetrahedral head mesh by loading it into Gmsh. In case freeview is not available, the spm viewer will be opened to allow for a basic check of the results.

Further notes
--------------

* Mesh resolution can be controlled using the -v option, which allows setting the vertex density (nodes per mm²) of the surface meshes. As a standard, :code:`headreco` uses 0.5 nodes per mm², resulting in head meshes with around 4 million tetrahedra.
* After the head mesh creation, temporary files are deleted to save disk space. Adding :code:`--noclean` prevents this.
* Manual editing: Edit one or more of the binary masks stored in :file:`m2m_{subID}/mask_prep/`. Then run :code:`headreco surfacemesh subID` and :code:`headreco volumemesh subID` to re-create the head mesh based on the edited masks. Add :code:`--no-cat` to the surfacemesh step in case you did nott use CAT12. Note: When using CAT12, surfaces instead of voxel masks will be stored for GM and WM in the :file:`mask_prep/` folder. For now, these surfaces cannot be manually improved.
* Transformation from and to MNI space: Both positions and results such as the electric field can be transformed between MNI and subject space. Please see below for a description of the corresponding command line programs. The transformation is based on a non-linear whole-head registration of the T1 of the subject to the MNI template that is determined during the SPM12 segmentation procedure. The transformations are stored in the :file:`m2m_{subID}/toMNI/` directory. Subject space is defined by the qform set in the :file:`m2m_{subID}/{subID}_T1fs_conform.nii.gz`, which can be found in the same folder as the head mesh. 
* When something goes wrong, you can check the :file:`m2m_{subID}/headreco_log.html` file.
* For headreco to run, MATLAB needs to be configured such that typing :code:`matlab` on a terminal window will start it. This is already the case in most Windows installations. However, on Linux and MacOSX you might need to create a link to the :code:`matlab` executable somewhere in your system :code:`$PATH`. This can be done with the commands

* Linux

  .. code-block:: bash
  
     sudo ln -s /usr/local/MATLAB/R<VERSION>/bin/matlab /usr/local/bin/matlab
  
  \

* MacOSX

  .. code-block:: bash
  
     sudo ln -s /Applications/MATLAB_R<VERSION>.app/bin/matlab /usr/local/bin/matlab
  
  \

If MATLAB is not installed in the default location, you can find out where it is installed by typing in a MATLAB terminal

.. code-block:: matlab

   matlabroot


References
-----------

`Nielsen, J. D., Madsen, K. H., Puonti, O., Siebner, H. R., Bauer, C., Madsen, C. G., ..., and Thielscher, A. (2018). Automatic skull segmentation from MR images for realistic volume conductor models of the head: Assessment of the state-of-the-art. NeuroImage, 174, 587-598. <https://doi.org/10.1016/j.neuroimage.2018.03.001>`_

