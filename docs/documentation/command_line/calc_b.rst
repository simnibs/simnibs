.. _calc_B_doc:

calc_B
===============

Description
-------------

Calculates the magnetic field caused by TES stimulation. For MREIT/MRCDI

.. note:: When using this feature, please cite `Yazdanian, H., Saturnino, G. B., Thielscher, A., & Knudsen, K. (2020). Fast evaluation of the Biot-Savart integral using FFT for electrical conductivity imaging. Journal of Computational Physics, 109408. <https://doi.org/10.1016/j.jcp.2020.109408>`_

Usage example
--------------

1. Run a simulation with the :code:`'J'` field as output.

2. Run :code:`calc_B`

.. code-block:: bash

  calc_B -i path/to/simulation_with_J.msh -r m2m_{subID}/T1fs_conform.nii.gz -o B.nii.gz

\

The magnetic field will be sampled in the same space as the :file:`m2m_{subID}/T1fs_conform.nii.gz` and will be stored in :file:`B.nii.gz`.
