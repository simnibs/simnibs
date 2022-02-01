.. _head_modeling_tutorial:

Creating Head Models
=====================

.. note:: The following describes head model creation using *charm*. For head modeling with *headreco* or *mri2mesh* in older SimNIBS versions, please see :ref:`here <head_modeling_tutorial_old>`.

MRI Scan
---------

To create individualized models, SimNIBS **requires** a T1-weighted image. T2-weighted images are optional, but **recommended**.

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

3. Run the segmentation and meshing pipeline *charm*. For example, for the *ernie* dataset, run

  .. code-block:: bash
  
     charm ernie org/ernie_T1.nii.gz org/ernie_T2.nii.gz
  
  \
  The subject ID (subID) :code:`ernie` is given as first argument. Charm will create named a folder :file:`m2m_ernie` that contains the segmentation results and the head mesh for the FEM simulations. The input images are given as final arguments (first the T1, then the T2).

\

  Alternatively, the reconstruction can be run with only the T1w image as input, but this can result in a less accurate skull region:

  .. code-block:: bash
  
     charm ernie org/ernie_T1.nii.gz
  
  \
  
Please see :ref:`charm_docs` for details and available options.

4. Check the segmentation. Click on the final segmentation viewer in the results.html (to be found in the m2m-folder of the subject). The viewer shows the outlines of the reconstructed tissue compartments, enabling a visual check whether the outlines are accurate.

Troubleshooting
----------------

* Please see the tutorial :ref:`fixheadmodel_tutorial` in case manually fixes to the segmentation are needed.


Further Reading
---------------

For more information on head meshing, please see:


`Puonti O, Van Leemput K, Saturnino GB, Siebner HR, Madsen KH, Thielscher A. (2020). Accurate and robust whole-head segmentation from magnetic resonance images for individualized head modeling. Neuroimage, 219:117044. <https://doi.org/10.1016/j.neuroimage.2020.117044>`_
