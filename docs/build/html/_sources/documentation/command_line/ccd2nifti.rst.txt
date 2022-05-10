.. _ccd2nifti_doc:

ccd2nifti
==========

Description
------------

This command line utility function allows conversion of coil dipole definition files (.ccd) files to nifti-1 volumes describing the magnetic vector potential (A-field).
The conversion utilizes the fast multipole methods to efficiently calculate the A-field on a grid. The size and resolution of this grid is defined in the header of the ccd file.

Usage example
-------------

1. Expand a specific ccd-file: Open a terminal, go to the coil directory and run

.. code-block:: bash

  ccd2nifti -i Drakaki_BrainStim_2022/MagVenture_C-B70.ccd
This will expand the file MagVenture_C-B70.ccd into MagVenture_C-B70.nii.gz, if the file does not already exist, in which case it will skip the expansion unless you force overwrite with the -f flag.

2. Expand coil within the current path recursively: Open a terminal, go to the intended directory and run

.. code-block:: bash

  ccd2nifti -i .


Further notes
---------------

* Type :code:`ccd2nifti -h` for the command line help
