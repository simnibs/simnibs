.. _mri2mesh_docs:

mri2mesh
==========

.. warning:: *Mri2mesh* is deprecated. Please use :ref:`charm_docs` instead. 

Description
-------------

mri2mesh reconstructs a tetrahedral head mesh from T1- and T2-weighted structural MR images. It runs also with only a T1w image, but will create better skull segmentations when also a T2w image is available.

.. attention:: mri2mesh depends on `FreeSurfer <https://surfer.nmr.mgh.harvard.edu/>`_ (5.3.0 or newer) and `FSL <https://fsl.fmrib.ox.ac.uk/fsl/fslwiki>`_ (5.0.5 or newer), and is therefore *not* compatible with Windows. Please see :ref:`optional_deps` for instructions on how to set-up FreeSurfer and FSL

Usage example
--------------

1. Open a terminal and go to the directory of the “Ernie” example data set.
2. Run the reconstruction: 

  .. code-block:: text
  
     mri2mesh --all ernie org/ernie_T1.nii.gz org/ernie_T2.nii.gz
  
  \

  The argument :code:`--all` tells :code:`mri2mesh` to run all reconstruction steps including volume meshing. The subject ID (subID) **ernie** is given next. Mri2mesh will create a mesh named :file:`ernie.msh`, a folder :file:`fs_ernie/` that contains the FreeSurfer results, and a folder :file:`m2m_ernie/` that contains the files that are needed for volume meshing. The input images are given as final arguments (first the T1, then the T2). When calling :code:`mri2mesh --all` the first time for a dataset, it will run FreeSurfer on it using the T1 as input. This is quite time-consuming. When re-running :code:`mri2mesh --all` it will use the existing FreeSurfer results, shortening the time required to ~3-4 hours.

  Alternatively, the reconstruction can be run with only the T1w image as input, but this will result in a less accurate skull region:

  .. code-block:: text
  
    mri2mesh --all ernie org/ernie_T1fs.nii.gz
  
  \

3. Check the results:

 .. code-block:: text
 
   mri2mesh -c ernie
 
 \

  This will show the reconstructed surfaces overlaid over the MR images using freeview. The red lines indicate the final surfaces used for volume meshing, the yellow indicate the GM and WM surfaces created by FreeSurfer. A second freeview will show the subject T1 overlaid on the MNI template for a visual check of the registration accuracy. In addition, you should have a look at the tetrahedral head mesh by loading it into gmsh.

Further notes
--------------

* A quick check can be performed by looking at the final volume masks overlaid over the structural images in fslview: :code:`mri2mesh --qc ernie`
* As a standard, :code:`mri2mesh` uses 60000 triangles for each white matter surface, and the number of triangles for the other surfaces are scaled relative to this number. This results in a volume mesh of ~3.5 million tetrahedra. Alternatively, you can adjust the mesh resolution by setting :code:`--numvertices=<mynumber>`.
* After the head mesh creation, temporary files are deleted to save disk-space. Adding :code:`--nocleanup` prevents this.
* When setting :code:`--t2pial`, FreeSurfer will use the T2 image to improve the estimate of the pial surfaces (recommended only for high-res with T2 images 1mm iso voxel). 
* Manual editing: For improving the GM and WM surfaces after the first run of :code:`mri2mesh`, edit the FreeSurfer results as described on the FreeSurfer wiki. Then run :code:`mri2mesh` again with the :code:`--all` option, as stated above. :code:`mri2mesh` will use the edited FreeSurfer results to create a new head mesh. For improving the ventricles, cerebellum, csf, skull or skin surfaces, manually edit one or more of the binary masks stored in :file:`m2m_{subID}/mask_prep/`. Then run :code:`mri2mesh` again wit the :code:`--all` and :code:`--keep_masks` options. The latter option will prevent :code:`mri2mesh` from overwriting the edited masks.
* Transformation from and to MNI space: Both positions and results such as the electric field can be transformed between MNI and subject space. Please see below for a description of the corresponding command line programs. The transformation is based on a non-linear whole-head registration of the T1 of the subject to the MNI template, using FSL’s fnirt command. Even though fnirt was developed for registering the brain, usually acceptable results are achieved by :code:`mri2mesh` for the whole head. When a T1 without fat suppression is used as input, the bright skull might be warped into the skin or brain. Using the option :code:`–mnimaskskull` can prevent this. A skull mask will then be applied to down-weight the skull intensity. The transformations are stored in the *toMNI* subdirectory of the :file:`m2m_{subID}/` folder. Subject space is defined by the qform set in the :file:`{subID}_T1fs_conform.nii.gz`, which can be found in the same folder as the head mesh. 
* When something goes wrong, you can check the :file:`m2m_{subID}/mri2mesh_log.html` file.

References
-----------

`Windhoff, M., Opitz, A., and Thielscher, A. (2013). Electric field calculations in brain stimulation based on finite elements: an optimized processing pipeline for the generation and usage of accurate individual head models. Human brain mapping, 34(4), 923-935. <https://doi.org/10.1002/hbm.21479>`_

