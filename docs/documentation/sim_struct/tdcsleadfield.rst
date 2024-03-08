.. _tdcsleadfield_doc:

TDCSLEADFIELD
==============

Describes the simulations for a Leadfield. Required for optimization


Initialization
---------------

* **Python**

  .. code-block:: python

     from simnibs import sim_struct
     tdcs_lf = sim_struct.TDCSLEADFIELD()

  \

* **MATLAB**

  .. code-block:: matlab

     tdcs_lf = sim_struct('TDCSLEADFIELD');

  \ 


Attributes
-----------

* **fnamehead**: *string (Python)/character array (MATLAB)*

  * **Desctiption** Name of head mesh file (:file:`{subID}.msh` file)
  * **Example**: *Python/MATLAB*

  .. code-block:: matlab

     tdcs_lf.fnamehead = 'ernie.msh'

  \ 

* **pathfem**: *string (Python)/character array (MATLAB)*

  * **Desctiption** Name of output folder
  * **Example**: *Python/MATLAB*

  .. code-block:: matlab

     tdcs_lf.pathfem = 'tdcs_leadfield/'

  \ 


* **field**: *'E' or 'J', optional*

  * **Description**: Whether to calculate the electric field *'E'* or current
    density *'J'*.
  * **Default**: 'E'


* **eeg_cap**: *string (Python)/character array (MATLAB), optional*

  * **Description**: Name of *.csv* file with EEG electrode positions.
  * **Default**: Automatically finds the file :file:`subpath/eeg_positions/EEG10-10_UI_Jurak_2007.csv` based on **fnamehead** or **subpath**
  * **Note**: Only needs to be set by the user if not using the standard *.csv* cap file.

* **interpolation**: *'middle gm', None/[], or list/cell array of strings (Python/MATLAB), optional. Default: 'middle gm'*

  * **Description**: Where to interpolate fields

    * *'middle_gm'*: Interpolate fields in the middle Gray Matter surface obtained from the head segmentation. Default value
    * *None/[]*: Do not interpolate the field anywhere, just store it in the region defined by the **tissues** attribute.
    * *list/cell array of strings*: List of mesh files in *'.stl'*, *'.gii'*, *'.off'* or *'.msh'* format. The files will be loaded and the fields will be interpolated at the mesh nodes.

  * **Default**: *'middle_gm'*

  .. note:: Does not work for *headreco* models ran with the :code:`--no-cat` option.

\

* **tissues**: *list (python) or array (MATLAB), optional*

  * **Description**: Tissues numbers of where to record the electric field, in addition to **interpolate**. Mixing surfaces and volumes is not allowed.

  * **Default**: [1006] (i.e. eye surfaces)

  * **Example**: *Python*

  .. code-block:: python

     # Example: Record electric fields in gray and white matter
     tdcs_lf.tissues = [1, 2]

  \ 

* **interpolation_tissue**:  *list (python) or array (MATLAB), optional*

  * **Description**: Tissues numbers to use as a base for the interpolation procedure. Nodes in the interpolated surfaces outside of the region defined by **interpolation_tissue** will be extrapolated using nearest neighbor. Ignored if **interpolation** is set to *None/[]*. Must correspond to a volume.

  * **Default**: [2] (i.e. gray matter surface)

* **electrode**: *ELECTRODE structure, list/array of ELECTRODE structures, or None/'none' optional*

  * **Description**: Electrodes to be used. Typically small round electrodes. There are 3
    ways to set this variable:

    * :ref:`electrode_struct_doc` structure: Use the same electrode shape for each electrode
      defined in the cap
    * list of :ref:`electrode_struct_doc` structures: Each electrode in the cap file will have
      the shape of the corresponding entry in the list
    * list of :ref:`electrode_struct_doc` structures and **eeg_cap** set to *None* (Python only):
      will use the **centre** and **pos_ydir** attributes of the electrodes to place
      them. This allows to set up electrodes on your own, without using a eeg cap provided by SimNIBS.
    * None(Python) or 'none' (MATLAB): Use point electrodes located at the surface nodes closest to the electrode center

  * **Default**: Use 1 x 1cm round electrodes with 4mm thickness


* **cond**: *list/array of COND structures (Python/MATLAB), optional*
   
  :ref:`Follow this link <cond_attribute_doc>`.

* **anisotropy_type**: *'scalar', 'vn', 'dir' or 'mc', optional*

  :ref:`Follow this link <anisotropy_type_attribute_doc>`.

* **aniso_maxratio**: *float*

  :ref:`Follow this link <aniso_maxratio_doc>`.

* **aniso_maxcond**: *float*

  :ref:`Follow this link <aniso_maxcond_doc>`.

* **fname_tensor**:*string (Python)/character array (MATLAB), optional*

  * **Description**: Name of NIfTI file with conductivity tensors
  * **Default**: Automatically finds the file :file:`d2c_{subID}/dti_results_T1space/DTI_conf_tensor.nii.gz` based on **fnamehead**.
  * **Note**: Only needed for simulations with anisotropic conductivities. And only needs to be set by the user if a file other than the above is to be used.

* **solver_options**: *string (pytohn) / character array (MATLAB)*

  :ref:`Follow this link <solver_options_doc>`.

References
-------------

`Saturnino, G. B., Siebner, H. R., Thielscher, A., & Madsen, K. H. (2019). Accessibility of cortical regions to focal TES: Dependence on spatial position, safety, and practical constraints. NeuroImage, 203, 116183. <https://doi.org/10.1016/j.neuroimage.2019.116183>`_

