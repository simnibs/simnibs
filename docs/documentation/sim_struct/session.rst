.. _session_doc:

SESSION
========

This is the basis of a SimNIBS simulation.
It contains information about the head mesh, the simulations to be performed, the post-processing steps and the output folder.

Initialization
---------------

* **Python**

  .. code-block:: python

     from simnibs import sim_struct
     s = sim_struct.SESSION()

  \

* **MATLAB**

  .. code-block:: matlab

     s = sim_struct('SESSION');

  \

Attributes
----------

* **fnamehead**: *string (Python)/character array (MATLAB)*

  * **Description** Name of head mesh file (given that you are in the `m2m_ernie` folder then :file:`{subID}.msh` file)
  * **Example**: *Python/MATLAB*

  .. code-block:: matlab

     s.fnamehead = 'ernie.msh'

  \

* **pathfem**: *string (Python)/character array (MATLAB)*

  * **Description** Name of output folder
  * **Example**: *Python/MATLAB*

  .. code-block:: matlab

     s.pathfem = 'simulation_outputs/'

  \


* **fields**: *string (Python)/character array (MATLAB), optional*

  * **Description**: Fields to be output. Any combination of

    * 'v' for voltage (in V)
    * 'e'/'E' for electric field magnitude/vector (in V/m)
    * 'j'/'J' for current density magnitude/vector (int A/m2).

  * **Default**: 'eE'
  * **Example**: *Python/Matlab*

  .. code-block:: matlab

     s.fields = 'eEjJ'

  \

.. _session_poslist:

* **poslist**: *list (Python)/cell array (MATLAB)*

  * **Description**: List or either :ref:`TDCSLIST <tdcslist_doc>` or :ref:`TMSLIST <tmslist_doc>` data structures
  * **Note**: In Python, this attribute does not need to be accessed directly. The user can use the *add_tdcslist* and *add_tmslist* methods instead.

* **open_in_gmsh**: *bool, optional*

  * **Description**: Whether to open the simulation results in *Gmsh* after finishing the simulation
  * **Default**: True (Python)/ False (MATLAB)
  * **Reference**: :ref:`Visualization tutorial <visualization_tutorial>`


* **map_to_surf**: *bool, optional*

  * **Description**: Whether to map the fields to the middle gray matter surface.
  * **Default**: False
  * **Reference**: `SimNIBS 2.1 tutorial paper <https://doi.org/10.1101/500314>`_

  .. warning:: Does not work in :ref:`headreco_docs` models ran with the :code:`--no-cat` option.

\

* **map_to_fsavg**: *bool, optional*

  * **Description**: Whether to map the fields to the FsAverage template.
  * **Default**: False
  * **Reference**: `SimNIBS 2.1 tutorial paper <https://doi.org/10.1101/500314>`_

  .. warning:: Does not work in :ref:`headreco_docs` models ran with the :code:`--no-cat` option.

\

* **map_to_vol**: *bool, optional*

  * **Description**: Whether to map the fields to a NIfTI volume. The NifTI volume will be in the same space as the :file:`m2m_{subID}/T1fs_conform.nii.gz` file.
  * **Default**: False
  * **Reference**: `SimNIBS 2.1 tutorial paper <https://doi.org/10.1101/500314>`_

* **map_to_mni**: *bool, optional*

  * **Description**: Whether to map the fields to the MNI template using a non-linear transformation.
  * **Default**: False
  * **Reference**: `SimNIBS 2.1 tutorial paper <https://doi.org/10.1101/500314>`_

* **tissues_in_niftis**: *list of tissue indices, or str ('all'), optional*

  * **Description**: List of the tissues for which the data will be interpolated to the niftis. To get fields everywhere, set to 'all'.
  * **Default**: 2 (i.e. gray matter)

* **subpath**: *string (Python)/character array (MATLAB), optional*

  * **Description**: Name of :file:`m2m_{subID}*` folder
  * **Default**: Automatically finds the :file:`m2m_{subID}/` folder based on **fnamehead**.
  * **Note**: Only required when **map_to_surf**, **map_to_fsavg**, **map_to_vol** or **map_to_mni** are set to *true*. Only needs to be set by the user if the *.msh* file was moved or renamed

* **fname_tensor**:*string (Python)/character array (MATLAB), optional*

  * **Description**: Name of NIfTI file with conductivity tensors
  * **Default**: Automatically finds the file :file:`d2c_{subID}/dti_results_T1space/DTI_conf_tensor.nii.gz` based on **fnamehead**.
  * **Note**: Only needed for simulations with anisotropic conductivities. And only needs to be set by the user if a file other than the above is to be used.

* **eeg_cap**: *string (Python)/character array (MATLAB), optional*

  * **Description**: Name of *.csv* file with EEG electrode positions
  * **Default**: Automatically finds the file :file:`subpath/eeg_positions/EEG10-10_UI_Jurak_2007.csv` based on **fnamehead** or **subpath**
  * **Note**: Only needed when setting TES electrodes of TMS coil positions using EEG electrode names. Only needs to be set by the user if not using the standard *.csv* cap file.



