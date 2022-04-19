.. _mni2subject_docs: 

mni2subject
===========

Description
------------

Transforms volumes from MNI space to subject space.

Usage example
---------------

1. Open a terminal and go to the directory of the “Ernie” example data set.
2. Suppose you have a volume mask defined in MNI space, called *Mask_MNI.nii.gz*
3. Run

.. code-block:: bash

  mni2subject -i Mask_MNI.nii.gz -m m2m_ernie/ -o Mask_ernie.nii.gz

\

4. This will create the *Mask_ernie.nii.gz* file, with the mask transformed to the subject space

Further notes
----------------

* Type :code:`mni2subject -h` for more information and options



