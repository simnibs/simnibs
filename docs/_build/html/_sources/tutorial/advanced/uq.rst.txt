.. _uq_tutorial:

Uncertainty Quantification
==========================


Tissue conductivity values are uncertain. Given that, it might be desirable to take this variability into consideration when simulating electric fields.

For that reason, we introduced a way in SimNIBS to easily perform Uncertainty Quantification (UQ) via Generalized Polynomial Chaos Expansion (gPC) for evaluating the effect of uncertainty in conductivity values on the fields of interest.

Before proceeding through this tutorial, please see `Saturnino et al., 2018 <https://doi.org/10.1016/j.neuroimage.2018.12.053>`_.


Introduction
-------------

In UQ via gPC, we quantify the uncertainty of the input variables (the conductivities) using a probability distribution.
SimNIBS supports two types of conductivity distribution: Normal and Beta.

After defining the uncertainty in the input variables, SimNIBS performs the UQ via gPC.
In gPC, we build a polynomial representation of the output variable (such as the electric field), given the input variable.
For SimNIBS, we developed an adaptive approach to select the polynomials to obtain fast convergence.
At each step, the error is evaluated using a K-means cross-validation scheme. The iteration stops when a tolerance is reached.



For more information on methodology, please see the supplementary material in  `Saturnino et al., 2018 <https://doi.org/10.1016/j.neuroimage.2018.12.053>`_.


Setting-up a Simulation with UQ
--------------------------------


It is simple to set-up a simulation with UQ. All you need to do is to set-up the  **distribution_type** and **distribution_parameters** parameters in the :ref:`cond_struct_doc` structures.

.. note:: As of now, UQ is only supported through scripting. For more information on scripts, please see :ref:`scripting_tutorial`.

.. warning:: Because conductivities can not be negative, we highly recommend the usage of beta-distributed random variables to represent conductivity uncertainties.

Python
''''''

.. literalinclude:: ../../../simnibs/examples/uncertainty_quantification/uncertainty_quantification.py
   :language: python


MATLAB
''''''

.. literalinclude:: ../../../simnibs/examples/uncertainty_quantification/uncertainty_quantification.m
   :language: matlab




Output Files
-------------
In the output folder of the examples above (the :file:`tdcs_uq/` folder) we have the output files

* :file:`ernie_TDCS_1_electrodes.msh`
    Head mesh (*Gmsh* format) with the electrodes
* :file:`ernie_TDCS_1_gpc.msh`
    Uncertainty quantification result in *Gmsh* format. There we have the data sets:

  * *{Field Name}_mean*
      Mean of the probability distribution describing the uncertainty of the given field
  * *{Field Name}_std*
      Standard Deviation of the probability distribution describing the uncertainty of the given field
  * *{Field Name}_sensitivity_{N}*
      Derivative-based sensitivity for the tissue number *N*. See `Saturnino et al., 2018 <https://doi.org/10.1016/j.neuroimage.2018.12.053>`_
  * *{Field Name}_sobol_{N}*
      Sobol coefficient for the tissue number *N*. See `Saturnino et al., 2018 <https://doi.org/10.1016/j.neuroimage.2018.12.053>`_

* :file:`ernie_TDCS_1_gpc.hdf5`: *HDF5* file containing partial information on the gPC expansion and partial results. The main datasets there are:

  * :file:`gpc_object`
      Information to reconstruct the gPC regressions
  * :file:`cond`
      Conductivity information
  * :file:`mesh`
      Original mesh (or with electrodes for tDCS simulations). Has the :file:`mesh/data_matrices/v_samples` data set with the electric potential values at all mesh nodes for all iterations of the adaptive algorithm.
  * :file:`mesh_roi`
      Mesh cropped to the ROI. Has the :file:`mesh/data_matrices/{Field Name}_samples` data sets with the field values (at nodes for *v* or elements for other quantities) for all iterations of the adaptive algorithm.

  .. note:: meshes stored in HDF5 files can be read with the :meth:`simnibs.msh.Msh.read_hdf5` class or with the *mesh_load_hdf5* function in MATLAB 


More Options
-------------

In *Python*, it is also possible to call lower-level functions to set more options for the UQ run and the post-processing of results.

In the example below, we set-up a UQ TMS problem with the ROI being the whole brain (tissues 1 and 2) and with a tolerance of 0.1.

.. literalinclude:: ../../../simnibs/examples/uncertainty_quantification/uq_setup_advanced.py
   :language: python


Secondary Quantities
---------------------

It is also possible to calculate secondary quantities, such as the 99th percentile of the electric field magnitude.

.. literalinclude:: ../../../simnibs/examples/uncertainty_quantification/uq_secondary_quantities.py
   :language: python

Acknowledgements
------------------

We would like to thank Konstantin Weise and Thomas Knoesche for the support in implementing the gPC in SimNIBS and supplying us with an early version of the `pygpc library <https://github.com/konstantinweise/pygpc>`_.

Further Reading
----------------
`Saturnino et al., A principled approach to conductivity uncertainty analysis in electric field calculations, Neuroimage 188, 2018 <https://doi.org/10.1016/j.neuroimage.2018.12.053>`_

`Weise, et al. Pygpc: A sensitivity and uncertainty analysis toolbox for Python, SoftwareX 11, 2020 <https://doi.org/10.1016/j.softx.2020.100450>`_
