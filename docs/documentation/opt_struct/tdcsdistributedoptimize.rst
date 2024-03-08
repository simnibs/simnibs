.. _tdcsdistributedoptimize_doc:


TDCSDistributedOptimize
=======================

TDCS optimization problem with a distributed target defined in a NIfTI file.

Initialization
---------------

* **Python**

  .. code-block:: python

     from simnibs import opt_struct
     opt = opt_struct.TDCSDistributedOptimize()

  \

* **MATLAB**

  .. code-block:: matlab

     opt = opt_struct('TDCSDistributedOptimize');

  \ 


Attributes
-----------

* **leadfield_hdf**: *string (Python)/character array (MATLAB)*

  * **Description**: Name of HDF5 file with leadfield (see :ref:`tdcsleadfield_doc` and also the first part of :ref:`tdcs_opt`)
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

* **target_image**: *str*

  * **Description**: Image to be emulated with the optimization
  * **Example**: *Python/MATLAB*

    Suppose we have a T-map representing a functional network called :file:`T_map.nii.gz`

  .. code-block:: matlab

     opt.target_image = 'T_map.nii.gz'

  \ 

* **mni_space**: *bool, optional*

  * **Descrption**: Wether the **target_image** is in MNI space.
  * **Default**: :code:`True`

* **subpath**: *str, optional*

  * **Descrption**: Path to the subject :file:`m2m_{SubID}` folder. Needed if :code:`mni_space=True`.
  * **Example**: *Python/MATLAB*

  .. code-block:: matlab

     opt.subpath = 'path/to/m2m_ernie'

  \ 

* **intensity**: *float, optional*
   
  * **Description**: Target field intensity. The larger the value, the more the optimization will hit the target network, but at the cost at having higher fields outside of it. Corresponds to :math:`E_0` in (`Ruffini et al., 2014 <https://doi.org/10.1016/j.neuroimage.2013.12.002>`_)
  * **Default**: 0.2

* **min_img_value**: *float ≥ 0, optional*

  * **Descrption**: minimum image value to be considered for optimization. All image values below it are set to zero. Corresponds to :math:`T_{min}` in (`Ruffini et al., 2014 <https://doi.org/10.1016/j.neuroimage.2013.12.002>`_)
  * **Default**: 0

* **open_in_gmsh**: *bool, optional*

  * **Descrption**: Whether to open the result in Gmsh after the calculations.
  * **Default**: :code:`True`

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



