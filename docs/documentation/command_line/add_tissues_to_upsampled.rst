.. _add_tissues_to_upsampled_doc:

add_tissues_to_upsampled
===========================

Description
------------

Adds extra tissues from an existing tissue label file to tissue_labeling_upsampled.nii.gz. The tissues to-be-added should be in a volume file, which is in the same space as the T1w scan used to create the segmentation.

Usage example
-------------
1. Download an example extra tissue file :download:`here <../../data/simnibs.nii.gz>`
2. Move the nifti file into the m2m-folder of the "Ernie" example data set.
3. Open a terminal and go to the directory of the “Ernie” example data set.
4. Make a copy of the tissue_labeling_upsampled.nii.gz, e.g., tissue_labeling_upsampled_orig.nii.gz
5. Run

.. code-block:: bash

  add_tissues_to_upsampled -i simnibs.nii.gz -t ./label_prep/tissue_labeling_upsampled.nii.gz -o tissue_labeling_upsampled.nii.gz --offset 50


6. The new *tissue_labeling_uspsampled.nii.gz* file will include the extra labels starting from label number 51. The offset values is added to the label values inside the extra tissue volume.

Further notes
---------------

* Type :code:`add_tissues_to_upsampled -h` for more information and options
* Please see :ref:`meshing_tutorial` for further information



