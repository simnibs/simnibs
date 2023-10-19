.. _analysis_tutorial:

Analyzing Simulations
======================

SimNIBS offers tools to help analyzing your simulation results.
Most of the tools have both *Python* and *MATLAB* versions. However, their arguments and outputs might be slightly different.

To run the analysis scripts, please follow the same steps needed to run the :ref:`simulation scripts <run_scripts>`.

.. _roi_mni_analysis:

ROI Analysis Using MNI Coordinates
-----------------------------------

The objective of this example is to evaluate the mean electric field in an ROI, defined using MNI coordinates.


This analysis follows the steps:

1. Load a simulation output.

        The default SimNIBS output is a **mesh**, that is, many nodes scattered in space forming tetrahedra and triangles. Electric fields (and also current density fields) are defined in **each tetrahedra and triangle**.

2. Select the gray matter volume elements for analysis:
        
        * *Python*:
                 :meth:`simnibs.Msh.crop_mesh` method
        * *MATLAB*:
                **mesh_extract_regions** function


3. Define the ROI.

   3.1. Use the **mni2subject_coords** function to transform coordinates from MNI to subject space.

   3.2. Get the center of the tetrahedra:

        * *Python*:
                :meth:`simnibs.Msh.elements_baricenters` method

        * *MATLAB*:
                **mesh_get_tetrahedron_centers** function

   3.3. Create a mask to separate elements inside the ROI.


4. Calculate field in ROI using a weighted average.

   4.1. Get the volumes of the tetrahedra:

        * *Python*:
                :meth:`simnibs.Msh.elements_volumes_and_areas` method


        * *MATLAB*:
                **mesh_get_tetrahedron_sizes** function

   4.2. Get the field of interest:

        * *Python*:
                :meth:`simnibs.Msh.field` attribute


        * *MATLAB*:
                **get_field_idx** function



   4.3. Calculate a weighted average:

Python
^^^^^^^

.. literalinclude:: ../../simnibs/examples/analysis/roi_analysis_mni.py
   :language: python


MATLAB
^^^^^^^
.. literalinclude:: ../../simnibs/examples/analysis/roi_analysis_mni.m
   :language: matlab


ROI Analysis Using Surfaces
---------------------------


An alternative analysis method is to use the **middle cortical surfaces** created by SimNIBS.


.. note:: This analysis method requires that the simulation has been interpolated to a cortical surface.
          This can be done either by setting the *Interpolate to cortical surface* option in the :ref:`sim_opt`. (under *Edit* -> *Simulation Options*),
          or by setting the **map_to_surf** option in the :ref:`session_doc` structure.

For this kind of analysis, we have to

1. Load the :file:`subject_overlays/*_central.msh` file.

   * This file contains both *left* and *right* hemispheres joined together in that order, as well as all the output fields.

   * Differently to the original mesh, this file has the fields defined in the **nodes**

   * SimNIBS also outputs each hemisphere as a *FreeSurfer surface* file and individual fields as a *FreeSurfer curv* file.

2. Define the ROI.

   2.1. Load an atlas transformed to subject space:

       * *Python*:
           :meth:`simnibs.subject_atlas` function
       * *MATLAB*:
           **subject_atlas** function

   2.2. Select a region from the atlas.

3. Calculate field in ROI using a weighted average.

   3.1. Calculate node areas:

      * *Python*:
           :meth:`simnibs.Msh.nodes_areas` method
      * *MATLAB*:
           **mesh_get_node_areas** function

   3.2. Get the field of interest:

       * *Python*:
           :meth:`simnibs.Msh.field` attribute


       * *MATLAB*:
           **get_field_idx** function



   3.3. Calculate a weighted average.

Python
^^^^^^^
.. literalinclude:: ../../simnibs/examples/analysis/roi_analysis_surf.py
   :language: python


MATLAB
^^^^^^^
.. literalinclude:: ../../simnibs/examples/analysis/roi_analysis_surf.m
   :language: matlab

