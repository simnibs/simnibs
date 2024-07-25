.. _tms_flex_opt:


General-purpose TMS Optimization
================================

.. note:: When using this feature in a publication, please cite 

\

Introduction
--------------
For setting up a general-purpose, flexible TMS optimization, you need to select the head model and the coil. In addition, you have to indicate whether your aim is to place a coil as close as possible to a target position on the head surface, or whether you want to optimize its placement to maximize the electric field magnitude in a region-of-interest (ROI) in the brain. For the first case, you will have to provide a target position *on the head surface* to which the coil center will be put as close as possible. For the second case, you need to provide a target region-of-interest (ROI) *in the brain*. In both cases, the coil position (and shape in case of flexible coils) will be chosen to prevent intersections of the coil casing with the head.

Example 1: Optimizing the electric field magnitude in a ROI
-----------------------------------------------------------
The following example optimizes the postion of a figure-8 coil to maximize the electric magnitude in a ROI around the left hand knob. As standard, a distance between the coil casing and head of 4 mm is ensured (as rough estimate of the distance caused by the hair). Here, it is set to 0 mm as example.

* *Python*

  .. literalinclude:: ../../simnibs/examples/tms_flex_optimization/tms_flex_fig8coil_emag.py
     :language: python

\

* *MATLAB*

  .. literalinclude:: ../../simnibs/examples/tms_flex_optimization/tms_flex_fig8coil_emag.m
     :language: matlab

\

The optimization will create the Gmsh output file :file:`ernie_MagVenture_MCF-B65_emag-optimization_surface_mesh.msh` with the ROI, the optimized field and the coil position:

.. figure:: ../images/Fig8_opt.png
   :scale: 40 %


Example 2: Optimizing the position and shape of a flexible coil
---------------------------------------------------------------
The following example optimizes the postion of a Brainsway H1 coil so that its center (as defined by the company) is as close as possible to a scalp position above the left DLPFC. In this case, the target coil position is *not* orthogonal to the local skin orientation underneath the coil center (i.e. the coil position cannot be defined by a target position and a second position indicating the y-direction). Therefore, the coil position has been here defined as 4x4 matrix in MNI space.

* *Python*

  .. literalinclude:: ../../simnibs/examples/tms_flex_optimization/tms_flex_Brainsway_H1_distance.py
     :language: python

\

* *MATLAB*

  .. literalinclude:: ../../simnibs/examples/tms_flex_optimization/tms_flex_Brainsway_H1_distance.m
     :language: matlab

\

The optimization will create the Gmsh output file :file:`ernie_Brainsway_H1_distance-optimization_head_mesh.msh` with the optimized field and coil position:

.. figure:: ../images/H1_opt.png
   :scale: 40 %
   

Notes
--------------
* Setting a suited starting position is recommended for flexible coils also when maximizing the electric field magnitude in a ROI.
* When maximizing the electric field magnitude in a ROI, the general-purpose TMS optimization uses the MKL Pardiso direct solver for accelerating the simulations. The SimNIBS standard FEM solver can be chosen optionally to reduce memory consumption, but will also substantially slow down the optimization.
* 32GB main memory are recommended, even thougth some optimizations will run with 16GB main memory.
* A combination of global and local search with settings that balance efficiency with robustness in finding a good solution is used as standard. For non-flexible coils, disabling global search (setting parameter run_global_optimization to False) will work fine for most situations. In case a more exhaustive optimization is desired, we suggest to set the "locally_biased" argument of the DIRECT solver to False.
* Please see :ref:`tms_flex_opt_doc` for a description of the option settings.
