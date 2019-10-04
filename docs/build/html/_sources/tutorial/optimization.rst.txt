TES Optimization
=================

SimNIBS 3 is able of optimizing TES montages for maximal focality.
There are two steps for performing the optimizations.

1. Create a *leadfield*. 
2. Set-up an optimization based on the *leadfield*.


Leadfield Calculations
-----------------------

The leadfield is calculated by first placing many electrodes in the head accordingly with an EEG cap, for example.
Afterwards, we calculate the field caused by activating each electrode individually, keeping a constant return electrode.
With this simulations, we form a matrix which can then be used to calculate the electric field caused by any combination of electrodes.
As the leadfield involves calculating many electric fields, it takes some time to run.

To calculate leadfields, the user must write a small Python of Matlab script using the :ref:`tdcsleadfield_doc` structure.

* *Python*

  .. literalinclude:: ../../simnibs/examples/calculate_leadfield.py
     :language: python

\

* *MATLAB*

  .. literalinclude:: ../../simnibs/examples/calculate_leadfield.m
     :language: matlab

\


Here, we are calculating a leadfield in the :file:`ernie` example data set, based on the EEG 10-10 system.
We use the :code:`map_to_surf` option to interpolate the electric field to the middle gray matter layer as this facilitates the optimization step.

.. note:: For more options in running tDCS leadfields, please see the documentation at :ref:`tdcsleadfield_doc`.

\

.. note:: the :code:`map_to_surf` option does not work with *headreco* models ran with the :code:`--no-cat` option.

\

This script will create the files in the :file:`tdcs_leadfield` folder:

* :file:`ernie_electrodes_EEG10-10_UI_Jurak_2007.msh`: File with the head mesh and electrodes placed on top. Can be used for checking electrode positions. Beware that some electrodes will be of the same color as the skin.

* :file:`ernie_ROI.msh`: The region of interest where the electric field values were recorded. In the example, the middle gray matter surface.

* :file:`ernie_leadfield_EEG10-10_UI_Jurak_2007.hdf5`: HDF5 file with the leadfield matrix, orgenized with the datasets:

  * :file:`mesh_electrodes`: Full mesh with the electrodes
  
  * :file:`mesh_leadfield`: Mesh with the region of interest only (in the example, the middle gray matter surface). Has the :file:`mesh_leadfield/leadfields/tdcs_leadfield` data set with the leadfield. In the attributes, you can se information such as the reference electrode, the field (*E* or *J*), units, currents and so on.

  .. note:: meshes stored in HDF5 files can be red with the :meth:`simnibs.msh.Msh.read_hdf5` class of with the *mesh_load_hdf5* function in MATLAB 


Optimization
-------------

Now, we will use the leadfield to optimize the electric field at a given target.

Safety constraints
~~~~~~~~~~~~~~~~~~~

The first step is to set :ref:`tdcsoptimize_doc` structure.
In this structure, we need to select the leadfield we will use for the optimization, a name for the optimization problem, safety constraints and limit the number of electrodes.

* *Python*

    .. code-block:: python
    
         from simnibs import optimization
         # Initialize structure
         opt = optimization.TDCSoptimize()
         # Select the leadfield file
         opt.leadfield_hdf = 'tdcs_leadfield/ernie_leadfield_EEG10-10_UI_Jurak_2007.hdf5'
         # Select a name for the optimization
         opt.name = 'tdcs_leadfield/optimization_example'
         # Select a maximum total current (in A)
         opt.max_total_current = 2e-3
         # Select a maximum current at each electrodes (in A)
         opt.max_individual_current = 1e-3
         # Select a maximum number of active electrodes (optional)
         opt.max_active_electrodes = 8


\

* *MATLAB*

    .. code-block:: matlab
    
         % Initialize structure
         opt = opt_struct('TDCSoptimize');
         % Select the leadfield file
         opt.leadfield_hdf = 'tdcs_leadfield/ernie_leadfield_EEG10-10_UI_Jurak_2007.hdf5';
         % Select a name for the optimization
         opt.name = 'tdcs_leadfield/optimization_example';
         % Select a maximum total current (in A)
         opt.max_total_current = 2e-3;
         % Select a maximum current at each electrodes (in A)
         opt.max_individual_current = 1e-3;
         % Select a maximum number of active electrodes (optional)
         opt.max_active_electrodes = 8;

\

.. note:: For more information see the :ref:`documentation for the TDCSoptimize strucure <tdcsoptimize_doc>`.

Target Settings
~~~~~~~~~~~~~~~

We now need to select the target for the optimization. The target positions are, as always in SimNIBS, given in **world coordinates** in **subject space** (:ref:`see here for more information <coords_doc>`). There are a few ways to obtain positions for optimizations in SimNBIS

* Open the :file:`{SubID}_ROI` file (eg.: :file:`tdcs_leadfield/ernie_ROI.msh`) in *Gmsh*. click in
        *Mesh* (on the left side of the window) -> *Inspect* and click in the model. (:download:`see here <../images/gmsh_inspect.png>`)

* Open the :file:`m2m_{subID}/T1fs_conform.nii.gz` file in a NifTi reader and

.. note:: The graphical user inteface is not an appropriate way to obtain positions for the optimization, as it only gets positions in the scalp surface, even when visualizing gray matter.

The input positions will be transformed by taking the closest point in the ROI.

Now, we write another script to run the optimization. We need target position and intensity

* *Python*

  .. code-block:: python

     target = opt.add_target()
     # Position of target
     target.positions = [-55.4, -20.7, 73.4]
     # Intensity of the electric field (in V/m)
     target.intensity = 0.2
     run_simnibs(opt)

\


* *MATLAB*

  .. code-block:: matlab

     % Position of target
     opt.target(1).positions = [-55.4, -20.7, 73.4];
     % Intensity of the electric field (in V/m)
     opt.target(1).intensity = 0.2;
     run_simnibs(opt)

\


As previously we ran the leadfield with the :code:`map_to_surf` setting, SimNIBS will optimize the electric field normal to the cortical surface. Otherwise, we should also set the **directions** attribute.

.. note:: For more options and information on targets please see the :ref:`referece for the TDCStarget structure <tdcstarget_doc>`.


Output files
'''''''''''''

The optimization outputs:

* :file:`{name}.csv`: comma separated values (CSV) files with optimal current values at each electrode (in A)
* :file:`{name}_electrodes.geo`: *Gmsh* *.geo* file for visualizing electrodes and currents
* :file:`{name}.msh`: *Gmsh* *.msh* file with the target and the optimized electric field in the ROI.
* :file:`{name}_summary.txt`: Some summary quantities about the optimization

Multiple targets
~~~~~~~~~~~~~~~~


To optimize multiple distant targets simultaneously, use multiple **target** structures.

* *Python*

  .. code-block:: python

     # Target in the left motor cortex
     target_left = opt.add_target()
     target_left.positions = [-55.4, -20.7, 73.4]
     target_left.intensity = 0.2
     # Target in the right motor cortex  
     target_right = opt.add_target()
     target_right.positions = [46.2, -35.8, 80.1]
     target_right.intensity = -0.2
     run_simnibs(opt)

\

* *MATLAB*

  .. code-block:: matlab

     % Target in the left motor cortex
     opt.target(1).positions = [-55.4, -20.7, 73.4];
     opt.target(1).intensity = 0.2;
     % Target in the right motor cortex  
     opt.target(2).positions = [46.2, -35.8, 80.1];
     opt.target(2).intensity = -0.2;
     run_simnibs(opt)

\

By using multiple targets, SimNIBS will try to hit each target with its intensity, whereas setting many **positions** in a single target, SimNIBS will try to hit the average intensity over the many positions.


Avoidance Regions
~~~~~~~~~~~~~~~~~~~


You can also add regions where the electric field should be more penalized. This is done using the **avoid** optional structure.

* *Python*

  .. code-block:: python

     avoid = opt.add_avoid()
     # Center of the region
     avoid.positions = [-35, -19, 85]
     # Radius of the region, in mm
     avoid.radius = 10

\


* *MATLAB*

  .. code-block:: matlab

     % Center of the region
     opt.avoid(1).positions = [-35, -19, 85];
     % Radius of the region, in mm
     opt.avoid(1).radius = 10;

\

.. note:: For more options and information on avoidance regions please see the :ref:`referece for the TDCSavoid structure <tdcsavoid_doc>`.



Complete example
-------------------


Click to download:

* :download:`Python <../../simnibs/examples/optimize.py>`
* :download:`MATLAB <../../simnibs/examples/optimize.m>`
