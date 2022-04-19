.. _optional_deps:

Optional Dependencies
======================

.. _matlab_setup:

MATLAB
-------

SimNIBS has a MATLAB API, available in the :file:`matlab_tools/` subfolder of the SimNIBS installation directory. Please add it to the matlab search path to make it available in your matlab.


FSL
----
**Required by:** :ref:`dwi2cond_docs`

1. Follow the instructions on `this link <http://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FslInstallation>`_

2. Add

.. code-block:: bash

   source /etc/fsl/5.0/fsl.sh

\
  in the end of the :file:`~/.bashrc` (Linux) or :file:`~/.bash_profile` (MacOSX) file (assuming that fsl is installed as usually into :file:`/etc/fsl/5.0`).
