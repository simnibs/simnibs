.. _tdcslist_doc:

TDCSLIST
========

Describes the electrode positions, shapes, currents as well as conductivity information for a single tDCS simulation.

Initialization
--------------

* **Python**

  .. code-block:: python

     from simnibs import sim_struct
     s = sim_struct.SESSION()
     tdcs_list = s.add_tdcslist()

  \

* **MATLAB**

  .. code-block:: matlab

     s = sim_struct('SESSION');
     tdcs_list = sim_struct('TDCSLIST');
     s.poslist{1} = tdcs_list;

  \ 

Attributes
-----------

* **currents** *list/array of floats (Python/MATLAB)*

  * **Description**: Current values (in **Ampere**), must sum to zero
  * **Example**: *Python/MATLAB*

  .. code-block:: python

     tdcs_list.currents = [1e-3, -1e-3]

  \

  * **Note**: Currents are given per channel, not per electrode. Please see the :ref:`ELECTRODE structure documentation <electrode_struct_doc>` for more information

* **electrode**: *list/array of ELECTRODE structures (Python/MATLAB)*

  * **Description**: List containing the electrodes to be used for the simulations.
  * **Note**: Please see the :ref:`ELECTRODE structure documentation <electrode_struct_doc>` for more information

.. _cond_attribute_doc:

* **cond**: *list/array of COND structures (Python/MATLAB), optional*

  * **Description**: List containing the tissue names and conductivities
  * **Default**: :ref:`Standard conductivity values <conductivity>`
  * **Example**:
    
    * *Python*

    .. code-block:: python

       # Change White Matter conductivity to 0.1 S/m
       # Notice, we need to reduce the tissue number by 1
       # in order to make up for the fact that python
       # indexing starts at zero
       tdcs_list.cond[0].value = 0.1

    \

    * *MATLAB*

    .. code-block:: matlab

       % Change White Matter conductivity to 0.1 S/m
       tdcs_list.cond(1).value = 0.1;

    \

  * **Note**: Please see the :ref:`COND structure documentation <cond_struct_doc>` for more information. All electrodes will get their conductivities from the tissues 100 and 500.

.. _anisotropy_type_attribute_doc:

* **anisotropy_type**: *'scalar', 'vn', 'dir' or 'mc', optional*

  * **Description**: Type of conductivity values to use in gray and white matter.

    * 'scalar': Isotropic, piecewise-constant conductivity values (default)
    * 'vn': Volume normalized anisotropic conductivities. In the volume normalization process, tensors are normalized to have the same trace and re-scaled according to their respective tissue conductivitiy (recommended for simulations with anisotropic conductivities, see `Opitz et al., 2011 <https://doi.org/10.1016/j.neuroimage.2011.06.069>`_)
    * 'dir': Direct anisotropic conductivity. Does not normalize individual tensors, but re-scales them accordingly to the mean gray and white matter conductivities (see `Opitz et al., 2011 <https://doi.org/10.1016/j.neuroimage.2011.06.069>`_).
    * 'mc': Isotropic, varying conductivities. Assigns to each voxel a conductivity value related to the volume of the tensors obtained from the direct approach (see `Opitz et al., 2011 <https://doi.org/10.1016/j.neuroimage.2011.06.069>`_).

  * **Default**: 'scalar'
  * **Note**: All options other than 'scalar' require conductivity tensors acquired from diffusion weighted images and processed with :ref:`dwi2cond <dwi2cond_docs>`.
  * **Reference**: `Tuch et al., 2001 <https://doi.org/10.1073/pnas.171473898>`_, `Opitz et al., 2011 <https://doi.org/10.1016/j.neuroimage.2011.06.069>`_

.. _aniso_maxratio_doc:

* **aniso_maxratio**: *float, optional*

  * **Description**: Maximum ratio between minimum an maximum eigenvalue in conductivity tensors
  * **Default**: 10
  * **Note**: Only taken into account when **anisotropy_type** is set to 'vn' and 'dir'
  * **Reference**: `Opitz et al., 2011 <https://doi.org/10.1016/j.neuroimage.2011.06.069>`_

.. _aniso_maxcond_doc:

* **aniso_maxcond**: *float, optional*

  * **Description**: Maximum mean conductivity value.
  * **Default**: 2 S/m
  * **Note**: Only taken into account when **anisotropy_type** is set to 'dir' or 'mc'
  * **Reference**: `Opitz et al., 2011 <https://doi.org/10.1016/j.neuroimage.2011.06.069>`_


.. _solver_options_doc:

* **solver_options**: *string (pytohn) / character array (MATLAB)*

  * **Description**: Options for the SimNIBS FEM solver.
  * **Default**: :code:`'-ksp_type cg -ksp_rtol 1e-10 -pc_type hypre -pc_hypre_type boomeramg -pc_hypre_boomeramg_coarsen_type HMIS'`
  * **Note**: Can be either a PETSc options string or the simple string :code:`'pardiso'` to use the MKL PARDISO solver.


