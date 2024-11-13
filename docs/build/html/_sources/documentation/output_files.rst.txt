.. _output_files:

Output Files
=============

Here, we will describe the output files of many SimNIBS see :ref:`file_formats` for a detailed description of each file type

Head Segmentation
------------------

:ref:`charm_docs` generates a folder :file:`m2m_<sub_id>/` with the segmentation results and the head mesh :file:`<sub_id>.msh` for the FEM simulations.


Simulations
-------------
The final simulation results are in the *Output Folder*, which can be set in the GUI or in scripts. By default, this folder is called :file:`simnibs_simulation/`, and is located in the same folder as the mesh file.

Below, we have:

* :file:`sub_id`: Subject ID (eg. ernie)
* :file:`poslist_nr`: Number of simulation in :ref:`poslist <session_poslist>`
* :file:`cond_type`: :ref:`Type of conductivity used <anisotropy_type_attribute_doc>`
* :file:`position_nr`: :ref:`Coil position number <tmslist_pos>`.
* :file:`coil_name`: Name of TMS coil model.

TDCS
''''

* :file:`<sub_id>_TDCS_<poslist_nr>_<cond_type>.msh`
    *Gmsh* mesh file with tDCS simulation results
* :file:`<sub_id>_TDCS_<poslist_nr>_el_currents.geo`
    *Gmsh* *.geo* file with simplified tDCS electrode models for visualization

TMS
'''

* :file:`<sub_id>_TMS_<poslist_nr>-<position_nr>_<coil_name>_<cond_type>.msh`
    *Gmsh* mesh file with TMS simulation results
* :file:`<sub_id>_TMS_<poslist_nr>-<position_nr>_<coil_name>_coil_pos.geo`
    *Gmsh* *.geo* file with a representation of the coil.

Others
'''''''

* :file:`simnibs_simulation_yyyymmdd.log`
    Simulation log
* :file:`simnibs_simulation_yyyymmdd.mat`
    Simulation configuration file, can be used to re-run the simulation
* :file:`fields_summary.txt`
    A few summary quantities for each field

Post-Processing
''''''''''''''''

In the section below, we will use 

* For tDCS simulations:
    :file:`<simulation_name>=TDCS_<poslist_nr>_<cond_type>`
* For TMS simulations:
    :file:`<simulation_name>=TMS_<poslist_nr>-<position_nr>_<coil_name>_<cond_type>` 

Surface Mapping
~~~~~~~~~~~~~~~~

* :file:`subject_overlays/`
    Surface overlays with fields in the subject specific space. Generated if the **Interpolate to cortical surface** (GUI) or the **map_to_surf** (script) are set


* :file:`fsavg_overlays/`
    Surface overlays with fields in FsAverage space. Generated if the **Interpolate to FsAverage surface** (GUI) or the **map_to_fsavg** (script) are set


Inside the folders, we have

* :file:`<sub_id>_<simulation_name>_<central/fsavg>.msh`
    Transformed result as a Gmsh *.msh* file for easy visualization
* :file:`<lh/rh>_<simulation_name>.<central/fsavg>.<field>.<quantity>`
    FreeSurfer generic overlay with given quantity (such as **magn** or **normal**) related to a given field.


Volume Mapping
~~~~~~~~~~~~~~

* :file:`subject_volumes/`
    NIfTI volumes with fields in a subject specific space. Generated if **Interpolate to a nifiti volume** (GUI) or **map_to_vol** (scripts) are set.

* :file:`mni_volumes/`
    NIfTI volumes with fields in MNI space. Generated if **Transform to MNI space** (GUI) or **map_to_mni** (scripts) are set.


Inside the folders, we have *NIfTI* files with the field values. For vector fields such as **E** and **J**, the *NIfTI* files are composed of 3 volumes, each corresponding to a direction
