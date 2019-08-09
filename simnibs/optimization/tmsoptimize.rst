.. _tmsoptimize_doc:

TMSoptimize
=============

Describes a TMS optimization problem

Initialization
---------------

* **Python**

  .. code-block:: python

     from simnibs import opt_struct
     opt = opt_struct.TMSoptimize()

  \

* **MATLAB**

  .. code-block:: matlab

     opt = opt_struct('TMSoptimize');

  \ 


Attributes
-----------

* **fnamehead**: *string (Python)/character array (MATLAB)*

  * **Desctiption**: Name of head mesh file (:file:`{subID}.msh` file)
  * **Example**: *Python/MATLAB*

  .. code-block:: matlab

     opt.fnamehead = 'ernie.msh'

  \ 

* **pathfem**: *string (Python)/character array (MATLAB)*

  * **Desctiption**: Name of output folder
  * **Example**: *Python/MATLAB*

  .. code-block:: matlab

     opt.pathfem = 'tms_optimization/'

  \ 

* **fnamecoil**: *string (Python)/ character array (MATLAB)*
  :ref:`Follow this link <tmslist_fnamecoil>`.

* **target**: *(3x1) list/array (Python/MATLAB) of floats*
  * **Desctiption**: Position of optimization target. Please see :ref:`coords_doc` for
    more information and how to obtain coordinates in SimNIBS
  * **Example**: *Python/MATLAB*

  .. code-block:: matlab

     opt.target = '[-43.4, -20.7, 83.4]'

  \ 

* **cond**: *list/array of COND structures (Python/MATLAB), optional*
   
  :ref:`Follow this link <cond_attribute_doc>`.

* **anisotropy_type**: *'scalar', 'vn', 'dir' or 'mc', optional*

  :ref:`Follow this link <anisotropy_type_attribute_doc>`.

* **aniso_maxratio**: *float, optional*

  :ref:`Follow this link <aniso_maxratio_doc>`.

* **aniso_maxcond**: *float, optional*

  :ref:`Follow this link <aniso_maxcond_doc>`.


* **target_size**: *float, optional*

  * **Description**: Radius of target area, in mm.
  * **Default**: 5
  
* **tissues**: *list/array (Python/MATLAB) of ints, optional*

  * **Description**: Tissues where the target is defined.
  * **Default**: [2] (Gray matter volume)

* **centre**: *(3x1) list/array (Python/MATLAB) of floats, optional*

  * **Description**: Position in scalp to use as a reference for the search space.
  * **Default**: *target* projected to the scalp

* **pos_ydir**: *(3x1) list/array (Python/MATLAB) of floats, optional*

  * **Description**: Reference position for the coil Y axis, with respect to the *pos* variable 
  * **Default**:  Search positions in a 360 degrees radius.

* **distance**: *float, optional*

  * **Description**: Distance from coil to scalp, in mm.
  * **Default**: 4

* **didt**: *float, optional*

  * **Description**: Coil dI/dt value, in A/s.
  * **Default**: 1e6

* **search_radius** *float, optional*

  * **Description**: Radius of search area, in mm.
  * **Default**: 20

* **spatial_resolution**: *float, optional*

  * **Description**: Space between coil positions, in mm.
  * **Default**: 5


* **angle_resolution**: *float, optional*

 * **Description**: Space between two coil angles, in degrees
 * **Default**: 30

* **search_angle**: *float, optional*

  * **Description**: Range of angles to use in search, in degrees.
  * **Default**: 360

* **open_in_gmsh**: *bool, optional*

  * **Description**: Wether to open the results in gmsh
  * **Default**: True
