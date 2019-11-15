.. _tdcstarget_doc:

TDCStarget
===========


Initialization
---------------

* **Python**

  .. code-block:: python

     from simnibs import optimization
     opt = optimization.TDCSoptimize()
     target = opt.add_target()



  \

* **MATLAB**

  .. code-block:: matlab

     opt = opt_struct('TDCSoptimize');
     opt.target(1)

  \ 


Attributes
-----------

.. _positions_attribute_doc:

* **positions**: *Nx3 list/array of floats (Python/MATLAB)*

  * **Desctiption**: Positions where the field is to be optimized. The positions
    are in **world coordinates** in **subject space** (:ref:`see here or more information about
    coordinates in SimNIBS <coords_doc>`). SimNIBS finds the position in the
    mesh closest to the input position. These can be obtained by

      * Open the :file:`tdcs_leadfield/{subID}_ROI.msh` file in *gmsh*, click in
        *Mesh* -> *Inspect* and click in the model
      * Open the :file:`m2m_{subID}/T1fs_conform.nii.gz` file in a NifTi reader and
        record the **world coordinates**.
      * Rea


* **directions**: *'normal' or Nx3 list/array of floats (Python/MATLAB), optional*

  * **Description**: Direction of the field to be optimized. If set to
    **normal** and the **leadfield** was run with **map_to_surf** set to *True*, will optimize the normal direction (see :ref:`tdcsleadfield_doc`).

  * **Defaut**: *'normal*

  .. note:: This argument is required if the leadfield obtained with :code:`map_to_surf = False`, which is the default in :ref:`tdcsleadfield_doc`


.. _indexes_attribute_doc:

* **indexes**: *Nx1 list/array of ints (Python/MATLAB), optional*

  * **Description**: As an alternative to **positions**, you can select the **node**
    index (if the leadfield was run with :code:`map_to_surf = True`, see :ref:`tdcsleadfield_doc`) or the **element** index
    (if the leadfield was run with :code:`map_to_surf = False`, see :ref:`tdcsleadfield_doc`)

  * **Default**: Get the points closest to the **positions**.


* **intensity**: *float, optional*

  * **Description**: Intensity of the field (*E* or *J*, see :ref:`tdcsleadfield_doc`) to
    be reached on average on the target and along the given direction. To optimize for
    intensity at the target rather than focality, you can set this value to a large
    number (eg: 100). With negative values, the direction will be inversed.
  * **Defaut**: 0.2


* **max_angle**: *float, optinal*

  * **Description**: Maximum angle between field and target direction.
  * **Default**: No maximum

  .. warning:: This condition is only fultilled in the mean, and not point-wise. Does not work in multi-target optimization


\

.. _radius_attribute_doc:

* **radius**: *float, optional*

 * **Description**: All points in the radius around the specified position/index to be added to the target area, in
   mm. Set to 0 if you want the target to be only the specified positions or indices.
 * **Default**: 2

* **tissues**: *list/array of ints (Python/MATLAB), optional*

  * **Descrption**: List of tissue indices where the target is defined. Leave empty if
    all tissues in the leadfield.
  * **Default**: All tissues

 
References
------------

`Saturnino, G. B., Siebner, H. R., Thielscher, A., & Madsen, K. H. (2019). Accessibility of cortical regions to focal TES: Dependence on spatial position, safety, and practical constraints. NeuroImage, 203, 116183. <https://doi.org/10.1016/j.neuroimage.2019.116183>`_
