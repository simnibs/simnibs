.. _gui_tutorial:

Setting up and Running Simulations
===================================

This tutorial is based on the :ref:`Example Dataset <dataset>`.
Please download it before continuing.

Starting the GUI and selecting a Head Model
--------------------------------------------
1. Launch the Graphical User Interface (GUI) from the Start Menu or by typing on a terminal window

  .. code-block:: bash
  
     simnibs_gui
  
  \

  A window such as the following should appear:

.. image:: ../images/tutorial_emptygui.png
   :align: center

2. Click on the *Browse* button next to the *m2m Folder* box.

.. image:: ../images/tutorial_bowse.png
   :align: center

3. Navigate to the example dataset folder and select the :file:`m2m_ernie` directory.

.. image:: ../images/tutorial_selectmsh.png
   :align: center

4. The head mesh will be loaded and the **m2m Folder** and **Output Folder** boxes filled.

.. image:: ../images/tutorial_guihead.png
   :align: center

\
  You can change the **Output Folder** by clicking on the *Browse* button close to it, and visualize the gray matter surface by clicking on the *Gray Matter* button above the head model. 


Setting up a tDCS Simulation
-----------------------------
1. Click on *Add tDCS Poslist* at the bottom left of the window.

.. image:: ../images/tutorial_addtdcs.png
   :align: center

\
  A new tab will show up.

2. Click on *Add Electrode*.

.. image:: ../images/tutorial_addelectrode.png
   :align: center

\
  A new row will be shown in the table.

3. Click twice in the *Shape* cell.

.. image:: ../images/tutorial_shapecell.png
   :align: center

\
  A new window will open where you can configure the electrode size and shape.

4. Let’s suppose we have a :math:`7.0 \times 5.5 \text{cm}^2` electrode, composed of a single gel layer of 5 mm thickness.

.. image:: ../images/tutorial_electrodegui.png
   :align: center

\
  Please refer to the :ref:`GUI documentation <gui_docs>` for a more detailed explanation of this window. You can also copy/paste electrode definitions by right-clicking in the shape cell.

5. Now double-click in the *Position* cell. A new window will appear. You can either select an EEG position from the drop-down menu in the left or click twice in the model to get a position.

.. image:: ../images/tutorial_position1.png
   :align: center

\

6. Choose the C3 electrode from the drop-down menu. A black sphere appears in the electrode position, with a green axis indicating the electrode's “y” axis.

.. image:: ../images/tutorial_position2.png
   :align: center

\


7. Select the electrode *current* in the first column of the electrode tables. **Positive** values designate **anodes**, and **negative** values **cathodes**. The sum of all current values has to be zero.

.. image:: ../images/tutorial_elcurrent.png
   :align: center

\

8. Now add a second electrode, place it at AF4 and double-click on a nearby position to rotate it, aligning the green axis with the top/down direction.

.. image:: ../images/tutorial_position3.png
   :align: center

\
  Select a shape and a current of :math:`-1.000 \text{mA}` for this second electrode.
  You can see a simple preview of the electrode shapes by clicking on *Preview Shapes* at the bottom of the screen.

Setting up a TMS Simulation
----------------------------
1. Click on *Add TMS Poslist* in the bottom of the window. A new TMS tab will be created.

.. image:: ../images/tutorial_addtms.png
   :align: center

\
2. Click on *Browse* and select the *Magstim_70mm_Fig8.ccd* coil file in the *legacy_and_other* subfolder

.. image:: ../images/tutorial_addcoil.png
   :align: center

\
3. Click on *Add Position*.

.. image:: ../images/tutorial_addposition.png
   :align: center

\
4. Double-click in the *Position* cell.

.. image:: ../images/tutorial_positioncell.png
   :align: center

\
5. You can select the coil position the same way as an electrode position, or alternatively by double-clicking on the head model. In the TMS case, the *y*-axis (green) indicates the direction of the coil handle. Here, we will click on the *Gray Matter* button above the model window and place the coil above the motor cortex, with the green axis pointing anteriorly. This means that the coil handle points posteriorly. :download:`See here for more information on coil coordinates <../data/coil_axesorientation.pdf>` .

.. image:: ../images/tutorial_coilpos.png
   :align: center

\

6. Additionally, you can also set the dI/dt (the current change ratio) and the coil-skin distance.

7. When using a *.nii.gz* coil file, click on *Show dA/dt field* to see the magnitude of the primary electric field.

.. image:: ../images/tutorial_dadt.png
   :align: center


\

.. attention:: This is **NOT the electric field**, but it can be interpreted as a very smooth approximation of it.

.. note:: Most coil files are supplied in *.tcd*/*.ccd*-format, which needs less disk space compared to *.nii.gz*. However, the preview option *Show dA/dt field* in the GUI currently works only for *.nii.gz* coil file. If needed, you can use the command line tool :ref:`coil2nifti_doc` to convert coil files from *.tcd*/*.ccd* to *.nii.gz*.

Setting Simulation Options
---------------------------
1. Go to *Edit* → *Simulation Options*.

.. image:: ../images/tutorial_simoptions.png
   :align: center


\
  The following window will appear:

.. image:: ../images/tutorial_simoptions2.png
   :align: center

\

2. We can select the *fields* to be output from the simulation.


  * **v**:
      Electrical Potential (Voltage). Units: Volts
  * **vector E**:
      Electric field vector. Units: V/m
  * **magn E**:
      Magnitude (or strength) of the electric field. Units: V/m
  * **vector J**:
      Current density vector. Units: A/m²
  * **magn J**:
      Magnitude of the current density. Units: A/m²
  * **Conductivities**:
      Conductivity field. For isotropic conductivities, this is a scalar.
      For anisotropic conductivities, this is the largest eigenvector of the conductivity tensor.
      Units: S/m
  * **dA/dt**:
      Primary field caused by the coil. TMS only. This is a vector field. Units: V/m

   Select **vector E** and **magn E**.

.. _tutorial_aditional_options:

3. We can also select *Additional Options*

  * **Open in Gmsh**:
      Opens the simulation results in *Gmsh*
  * **Interpolate to cortical surface**:
      Interpolates the fields along a surface at the center of the gray matter sheet. Not available for :ref:`headreco_docs` models ran with :code:`--no-cat`.
  * **Transform to fsaverage space**:
      Interpolates to the middle of gray matter and transforms it to FsAverage space. Not available for :ref:`headreco_docs` models ran with :code:`--no-cat`.
  * **Interpolate to a nifiti volume**:
      Interpolates the fields to a nifti volume.
  * **Transform to MNI space**:
      Interpolates the fields to a nifti volume and applies a transformation to MNI space.

  For the example run, we will select all of the above.


Running a Simulation
---------------------
1. Click on Run at the bottom of the screen.

.. image:: ../images/tutorial_runsim.png
   :align: center

\
2. If there are no errors in the problem set-up, a new window will appear and show the simulation progress. The simulation takes a few minutes, and when finished a Gmsh window opens with the simulation results.

Now, please go on to our tutorial on :ref:`visualization_tutorial`.

Output Files
-------------

After the simulation is finished, the :file:`simnibs_simulation` directory will look like the following:

.. image:: ../images/tutorial_output_files.png
   :align: center

\

The main files here are the :file:`.msh` files

  * :file:`ernie_TDCS_1_scalar.msh`
      With results of the tDCS simulation
  * :file:`ernie_TMS_2-0001_Magstim_70mm_Fig8_nii_scalar.msh`
      With results of the TMS simulation

The folders
   * :file:`fsavg_overlays`
   * :file:`mni_volumes`
   * :file:`subject_overlays`
   * :file:`subject_volumes`

are only present if the corresponding :ref:`options <tutorial_aditional_options>` are
selected.


For a complete explanation of the output, please see :ref:`output_files` .

Further Reading
----------------

For more information on the GUI, please see the `SimNIBS 2.1 tutorial paper <https://doi.org/10.1101/500314>`_.
