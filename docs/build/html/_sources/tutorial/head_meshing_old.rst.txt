.. _head_modeling_tutorial_old:

.. warning:: The following describes head model creation based on *headreco* or *mri2mesh*. Both are deprecated. Please use :ref:`charm_docs` instead. 

Creating Head Models
=====================

MRI Scan
---------

To create individualized models, SimNIBS **requires** a T1-weighted image. T2-weighted images are optional, but **highly recommended**.

T1-weighted images
~~~~~~~~~~~~~~~~~~~

Commonly used for segmentation of the brain. Usually acquired at a rather low readout bandwidth to maximize SNR. However, using a low readout bandwidth comes with the disadvantage that the positions of the (fatty) spongy bone and subcutaneous fat will be displaced in the MR images due to the chemical shift artifact. This can result in the spongy bone touching brain gray matter (GM), rendering an accurate segmentation of the GM pial surface and the boundary between cerebrospinal fluid (CSF) and skull difficult. For this reason, we recommend to **acquire the T1w images with a fat suppression method** such as selective water excitation. Specific sequence parameters can be seen at `Nielsen et al., NeuroImage, 2018 <https://doi.org/10.1016/j.neuroimage.2018.03.001>`_.


T2-weighted images
~~~~~~~~~~~~~~~~~~~
As both CSF and compact bone are dark in T1w images (also spongy bone will be dark when fat suppression is used), an accurate reconstruction of the skull is difficult when only a T1w image is available. For this reason, we recommend to additionally **acquire a T2w image without fat suppression**. T2w images are usually acquired at a high readout bandwidth, minimizing the chemical shift artifact. As they also provide very good contrast between CSF (bright) and compact bone (dark), they are a good starting point for segmenting the skull. For specific sequence parameters, please see `Nielsen et al., NeuroImage, 2018 <https://doi.org/10.1016/j.neuroimage.2018.03.001>`_.



Good and Bad T1w and T2w images
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. figure:: ../images/t1_t2.jpg

  **Good** T1w and T2w images  

\

 
.. figure:: ../images/bad_t1.png
   :scale: 50 %
                        
   **Bad** T1w image. The spongy bone is in part very bright (indicated by the red arrows) and a dark CSF region to separate GM and skull is missing at several positions (purple arrows indicate example positions).

\ 

Segmentation and Meshing
-------------------------

After scanning and having the MRI images in NifTI format, the next step is to create a head mesh. To do it,

1. Start a Command Prompt/Terminal window (:ref:`Windows Instructions <windows_terminal>`)

2. Navigate to the folder where the NifTI files are located (use the :code:`cd` command)

3. Run a segmentation and meshing pipeline. SimNIBS offers two tools for this job, :ref:`headreco_docs` and :ref:`mri2mesh_docs`. *headreco* is the most recent tool, and should be preferred in most cases.

.. list-table::
   :widths: 33 33 33
   :header-rows: 1

   * -
     - :ref:`headreco_docs`
     - :ref:`mri2mesh_docs`
   * - Operating Systems
     - Windows, Linux, MacOSX
     - Linux, MacOSX
   * - External Dependencies
     - MATLAB
     - FreeSurfer and FSL
   * - Coverage
     - Whole head and neck
     - Above the mouth
   * - Time to run
     - ~2 hours (with CAT12), ~1 hour (without CAT12)
     - ~10 hours
   * - Reference
     - `Nielsen et al., 2018 <https://doi.org/10.1016/j.neuroimage.2018.03.001>`_
     - `Windhoff et al., 2013 <https://doi.org/10.1002/hbm.21479>`_

\
  Before running, please read the documentation of the tool of your choice.

4. Check the segmentation. This can be done with :code:`headreco check <SUB_ID>` for models ran with *headreco* or :code:`mri2mesh -c <SUB_ID>` for models ran with *mri2mesh*. Please see the documentation for your tool of choice for more information. With *FreeSurfer* installed, this command opens two *freeview* windows, one with the T1w image, surface outlines and tissue masks, and another with the MNI transformed T1w image overlaid on the MNI template. For *headreco*, if FreeSurfer is not installed, a SPM window opens showing the MNI registration.

.. figure:: ../images/check_segmentation_good.png

   **Good** segmentation. Notice how the surface oulines (in white) nicely follows the tissue shapes in the T1w image

\

.. figure:: ../images/check_segmentation_bad.png

   **Bad** segmentation. Notice how part of the skull is missing

\

.. figure:: ../images/check_mni.png

   Example of a good MNI registration.

\

5. Load the head model in *Gmsh* and check the head mesh to ensure that head meshing went fine. First go to *Tools -> Options -> Mesh* and select *Volume faces*.  Then go to *Tools -> Visibility* and select the volumes one by one.

Troubleshooting
----------------

* Sometimes, Gmsh fails to mesh one on more tissue compartments. When this happens, running the whole pipeline again and increasing the mesh density with the :code:`-v` argument in :ref:`headreco_docs` or the :code:`--numvertices` argument in :ref:`mri2mesh_docs`.


Further Reading
---------------

For more information on head meshing, please see our `SimNIBS 2.1 tutorial paper <https://doi.org/10.1101/500314>`_, `Nielsen et al., 2018 <https://doi.org/10.1016/j.neuroimage.2018.03.001>`_ and `Windhoff et al., 2013 <https://doi.org/10.1002/hbm.21479>`_.


