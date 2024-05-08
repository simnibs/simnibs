.. _tdcsoptimize_doc:

TDCSoptimize
=============

Describes a tDCS optimization problem

Initialization
---------------

* **Python**

  .. code-block:: python

     from simnibs import opt_struct
     opt = optimization.TDCSoptimize()

  \

* **MATLAB**

  .. code-block:: matlab

     opt = opt_struct('TDCSoptimize');

  \ 


Attributes
-----------

.. _leadfield_hdf_doc:

* **leadfield_hdf**: *string (Python)/character array (MATLAB)*

  * **Desctiption**: Name of HDF5 file with leadfield (see :ref:`tdcsleadfield_doc`)
  * **Example**: *Python/MATLAB*

  .. code-block:: matlab

     opt.leadfield_hdf = 'tdcs_leadfield/ernie_leadfield_EEG10-10_UI_Jurak_2007.hdf5'

  \ 

.. _opt_name_doc:

* **name**: *string (Python)/character array (MATLAB)*

  * **Desctiption**: Name of the optimization problem. Gives the path and the prefix to
    the output files
  * **Example**: *Python/MATLAB*

    To have the output files from the optimization in the :file:`tdcs_leadfield/` folder
    and the :file:`motor_cortex` prefix

  .. code-block:: matlab

     opt.name = 'tdcs_leadfield/motor_cortex'

  \ 

* **target**: *list/array of TDCStarget structure (Python/MATLAB)*

  * **Description**: List of targets for optimization. See :ref:`tdcstarget_doc`. By
    setting a list of targets, SimNIBS will try to optimize all of them simultaneously.

* **avoid**: *list/array of TDCSavoid structure (Python/MATLAB), optional*

  * **Description**: List of regions to be more heavily punished during optimization. See :ref:`tdcsavoid_doc`.

* **max_total_current**: *float, optional*

  * **Description**: Maximum total injected current, in A
  * **Default**: 2e-3

* **max_individual_current**: *float, optional*

  * **Description**: Maximum current injected in each electrode, in A
  * **Default**: 1e-3

* **max_active_electrodes**: *int, optional*

  * **Description**: Maximum number of active electrodes. Leave empty if no maximum
    number of electrodes
  * **Default**: No maximum

References
-----------

`Saturnino, G. B., Siebner, H. R., Thielscher, A., & Madsen, K. H. (2019). Accessibility of cortical regions to focal TES: Dependence on spatial position, safety, and practical constraints. NeuroImage, 203, 116183. <https://doi.org/10.1016/j.neuroimage.2019.116183>`_


