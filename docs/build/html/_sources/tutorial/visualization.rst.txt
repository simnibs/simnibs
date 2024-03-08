.. _visualization_tutorial:

Visualizing Results
====================
This tutorial comes after the tutorial on :ref:`gui_tutorial`. Please be sure to complete it before going further.


Visualizing Results in Gmsh
----------------------------
`Gmsh <http://gmsh.info/>`_ is a powerful program for 3D visualization. However, it has a somewhat steep learning curve.
For this reason, SimNIBS writes *.opt* files with basic visualization configurations in addition to the mesh files.
Here we will give a walk-through of a few useful features in Gmsh.

.. note:: `Gmsh <http://gmsh.info/>`_ is distributed together with SimNIBS, you don't need to install it separately.

Starting Gmsh
~~~~~~~~~~~~~~

Gmsh automatically opens when a simulation is finished.
Alternatively, you can start Gmsh by double clicking in the :file:`ernie_TDCS_1_scalar.msh` file or typing in the terminal

  .. code-block:: bash
  
     gmsh ernie_TDCS_1_scalar.msh
  
  \
  The following window will appear:

.. image:: ../images/tutorial_gmshinit.png
   :align: center

\

  We see the gray matter surface with the **magnitude** of the electric field (**magnE**), as well as the electrode currents. You can rotate the model with the left mouse button, translate it with the right button, and zoom with the mouse wheel.

 
Selecting Where to Visualize
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
To visualize another surface, such as white matter:

1. Go to *Tools* → *Visibility*

.. image:: ../images/tutorial_gmshvisibility.png
   :align: center

\

 and select *Elementary entities*

.. image:: ../images/tutorial_elementaryentities.png
   :align: center

\

2. Select surfaces and volumes to be visualized. In SimNIBS, we have defined

   1: White Matter volume

   2: Gray Matter volume

   3: Cerebrospinal Fluid (CSF) volume

   4: Skull volume

   5: Skin volume

   6: Eye volumes (headreco models only)

   101+: Electrode rubber layer volume

   501+: Electrode gel / sponge volume

   The corresponding surfaces have the number 1000 + volume value. For example, 1002 corresponds to the gray matter surface. The only exception are the electrode contacts, numbered 2100+. You can select multiple surfaces or volumes by holding *Ctrl*. Selecting 1001 and clicking on *Apply*, we can see the white matter surface:

.. image:: ../images/tutorial_gmshwmfield.png
   :align: center

\


Selecting What to Visualize
~~~~~~~~~~~~~~~~~~~~~~~~~~~


You can select the field to be visualized by checking the boxes in the left of the Gmsh window.

.. image:: ../images/tutorial_gmshnorme.png
   :align: center

\

To visualize the surfaces or volumes without any field, go to *Tools* -> *Options* -> *Mesh* and click on *Surface Faces* or *Volume Faces*

.. image:: ../images/tutorial_gmshdeselecsurf.png
   :align: center

\

Changing the Scale
~~~~~~~~~~~~~~~~~~~
To change the scale of a field visualization, select *Tools* → *Options* → *View [N]* and change the *Min* and *Max* values.

.. image:: ../images/tutorial_gmshview.png
   :align: center

\


Exporting an Image
~~~~~~~~~~~~~~~~~~~

To create an image, go to *File* -> *Export* (or press *Ctrl+E*) and type in a file name with a  *.png* or  *.jpg* extension.


Other Functionalities
~~~~~~~~~~~~~~~~~~~~~~

We recommend users to explore the many functionalities of Gmsh. One can, for example, produce the image below by selecting the Volume 2 for visualization an clipping the model in *Tools* →  *Clipping*.

.. image:: ../images/tutorial_gmshclip.png
   :align: center

\


Visualizing Results in MATLAB
------------------------------

To visualize results in MATLAB, add the :file:`<SIMNIBS_INSTALL_DIR>/matlab/` folder to your MATLAB path. Afterwards, type

.. code-block:: matlab

   mesh_get_simulation_result

\


or use

.. code-block:: matlab

   mesh_show_surface

\

to have more control over the visualization.

Further Reading
----------------
For more information on visualization and simulation output, please see the `SimNIBS 2.1 tutorial paper <https://doi.org/10.1101/500314>`_.
