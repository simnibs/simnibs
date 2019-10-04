.. _tdcsoptimize_doc:

TDCSoptimize
=============

Describes a tDCS optimization problem

Initialization
---------------

* **Python**

  .. code-block:: python

     from simnibs import optimization
     opt = optimization.TDCSoptimize()

  \

* **MATLAB**

  .. code-block:: matlab

     opt = opt_struct('TDCSoptimize');

  \ 


Attributes
-----------

* **leadfield_hdf**: *string (Python)/character array (MATLAB)*

  * **Desctiption**: Name of HDF5 file with leadfield (see :ref:`tdcsleadfield_doc`)
  * **Example**: *Python/MATLAB*

  .. code-block:: matlab

     opt.leadfield_hdf = 'tdcs_leadfield/ernie_leadfield_EEG10-10_UI_Jurak_2007.hdf5'

  \ 

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

  * **Description**: Maximum total injected current, in mA
  * **Default**: 2e-3

* **max_individual_current**: *float, optional*

  * **Description**: Maximum current injected in each electrode, in mA
  * **Default**: 1e-3

* **max_active_electrodes**: *int, optional*

  * **Description**: Maximum number of active electrodes. Leave empty if no maximum
    number of electrodes
  * **Default**: No maximum


Examples
---------

* Run an optimization based on the :file:`tdcs_leadfield/ernie_leadfield_EEG10-10_UI_Jurak_2007.hdf5`.
  
  * *Python*

    .. code-block:: python
    
         from simnibs import optimization
         opt = optimization.TDCSoptimize()
         opt.leadfield_hdf = 'tdcs_leadfield/ernie_leadfield_EEG10-10_UI_Jurak_2007.hdf5'
         opt.name = 'tdcs_leadfield/motor_cortex'
         opt.max_total_current = 2e-3
         opt.max_individual_current = 1e-3

     \


