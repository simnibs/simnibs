.. _tdcs_opt:

TDCS Optimization
=================

SimNIBS 3 is able of optimizing TDCS montages for maximal focality.
There are two steps for performing the optimizations.

1. Create a *leadfield*. 
2. Set-up an optimization based on the *leadfield*.


Leadfield Calculations
-----------------------

The leadfield is calculated by first placing many electrodes in the head accordingly with an EEG cap, for example, and afterwards calculating the field caused by each electrode individually, keeping a constant return electrode.

With this simulations, we form a matrix which can then be used to calculate the electric field caused by any combination of electrodes.
As the leadfield involves calculating many electric fields, it takes some time to run, but usually no more than 1h.

To calculate leadfields, the user must write a small Python of Matlab script using the :ref:`tdcsleadfield_doc` structure. Unless specified differently, the leadfields are calculated for EEG10-10 electrode positions.

* *Python*

  .. literalinclude:: ../../simnibs/examples/optimization/leadfield.py
     :language: python

\

* *MATLAB*

  .. literalinclude:: ../../simnibs/examples/optimization/leadfield.m
     :language: matlab

\


Here, we are calculating a leadfield in the :file:`ernie` example data set, based on the EEG 10-10 system.

.. note:: For more options in running tDCS leadfields, please see the documentation at :ref:`tdcsleadfield_doc`.

\


This script will create the files in the :file:`leadfield` folder:

* :file:`ernie_electrodes_EEG10-10_UI_Jurak_2007.msh`: File with the head mesh and electrodes placed on top. Can be used for checking electrode positions. Beware that some electrodes will be of the same color as the skin.

* :file:`ernie_ROI.msh`: The region of interest where the electric field values were recorded. In the example, the middle gray matter surface.

* :file:`ernie_leadfield_EEG10-10_UI_Jurak_2007.hdf5`: HDF5 file with the leadfield matrix, orgenized with the datasets:

  * :file:`mesh_electrodes`: Full mesh with the electrodes
  
  * :file:`mesh_leadfield`: Mesh with the region of interest only (in the example, the middle gray matter surface). Has the :file:`mesh_leadfield/leadfields/tdcs_leadfield` data set with the leadfield. In the attributes, you can see information such as the reference electrode, the field (*E* or *J*), units, currents and so on.

  .. note:: meshes stored in HDF5 files can be read with the :meth:`simnibs.msh.Msh.read_hdf5` class or with the *mesh_load_hdf5* function in MATLAB 


Optimization
-------------

Now, we will use the leadfield to optimize the electric field at a given target.

Simple Optimiztion
~~~~~~~~~~~~~~~~~~~

The first step is to set :ref:`tdcsoptimize_doc` structure.
In this structure, we need to select the leadfield we will use for the optimization, a name for the optimization problem, safety constraints and limit the number of electrodes.

Afterwards, we need to define an optimization target using an :ref:`tdcstarget_doc` structure.

* *Python*

  .. literalinclude:: ../../simnibs/examples/optimization/tdcs_optimize.py
     :language: python

\

* *MATLAB*

  .. literalinclude:: ../../simnibs/examples/optimization/tdcs_optimize.m
     :language: matlab

\

.. note:: For more information see the documentation for :ref:`tdcsoptimize_doc` and :ref:`tdcstarget_doc`.

Output files
'''''''''''''

The optimization outputs:

* :file:`{name}.csv`: comma separated values (CSV) files with optimal current values at each electrode (in A)
* :file:`{name}_electrodes.geo`: *Gmsh* *.geo* file for visualizing electrodes and currents
* :file:`{name}.msh`: *Gmsh* *.msh* file with the target and the optimized electric field in the ROI.
* :file:`{name}_summary.txt`: Some summary quantities about the optimization


Maximizing intensity
~~~~~~~~~~~~~~~~~~~~~

To maximize intensity at the target, disregarding field focality, simply use a large value for the target intensity.


* *Python*

  .. literalinclude:: ../../simnibs/examples/optimization/tdcs_optimize_intensity.py
     :language: python

\

* *MATLAB*

  .. literalinclude:: ../../simnibs/examples/optimization/tdcs_optimize_intensity.m
     :language: matlab

\

Using MNI Coordinates 
~~~~~~~~~~~~~~~~~~~~~

The target positions are, as always in SimNIBS, given in **world coordinates** in **subject space** (:ref:`see here for more information <coords_doc>`). However, we can use the *mni2subject_coords* function to transform coordinates from MNI space to subject space. When the transformed coordinates are outside the gray matter of the subject, they will be projected to the closest gray matter position.


* *Python*

  .. literalinclude:: ../../simnibs/examples/optimization/tdcs_optimize_mni.py
     :language: python

\

* *MATLAB*

  .. literalinclude:: ../../simnibs/examples/optimization/tdcs_optimize_mni.m
     :language: matlab

\

Multiple targets
~~~~~~~~~~~~~~~~


To optimize multiple distant targets simultaneously, just use multiple **target** structures.

* *Python*

  .. literalinclude:: ../../simnibs/examples/optimization/tdcs_optimize_multi_target.py
     :language: python

\

* *MATLAB*

  .. literalinclude:: ../../simnibs/examples/optimization/tdcs_optimize_multi_target.m
     :language: matlab

\

By using multiple targets, SimNIBS will try to hit each target with its intensity, whereas setting many **positions** in a single target, SimNIBS will try to hit the average intensity over the many positions.


Avoidance Regions
~~~~~~~~~~~~~~~~~~~

You can also add regions where the electric field should be more penalized than elsewhere. This is done using the **avoid** optional structure.

* *Python*

  .. literalinclude:: ../../simnibs/examples/optimization/tdcs_optimize_avoid.py
     :language: python

\

* *MATLAB*

  .. literalinclude:: ../../simnibs/examples/optimization/tdcs_optimize_avoid.m
     :language: matlab

\


.. note:: For more options and information on avoidance regions please see the :ref:`referece for the TDCSavoid structure <tdcsavoid_doc>`. You can visualize the position of the avoided region in the results by deselecting **normE** in gmsh, and selecting **avoid_1**.


Avoid field in Eyes
~~~~~~~~~~~~~~~~~~~~

Here, we will perform an optimization penalizing more the electric field in the eyes

Calculate Leadfield with Eyes
''''''''''''''''''''''''''''''

First we need to calculate a Leadfield which contains the eyes.

* *Python*

  .. literalinclude:: ../../simnibs/examples/optimization/leadfield_eyes.py
     :language: python

\

* *MATLAB*

  .. literalinclude:: ../../simnibs/examples/optimization/leadfield_eyes.m
     :language: matlab

\


Optimize
'''''''''

Now we set-up the optimization and select the eyes as and avoid region

* *Python*

  .. literalinclude:: ../../simnibs/examples/optimization/tdcs_optimize_avoid_eyes.py
     :language: python

\

* *MATLAB*

  .. literalinclude:: ../../simnibs/examples/optimization/tdcs_optimize_avoid_eyes.m
     :language: matlab

\

References
------------

`Saturnino, G. B., Siebner, H. R., Thielscher, A., & Madsen, K. H. (2019). Accessibility of cortical regions to focal TES: Dependence on spatial position, safety, and practical constraints. NeuroImage, 203, 116183. <https://doi.org/10.1016/j.neuroimage.2019.116183>`_


