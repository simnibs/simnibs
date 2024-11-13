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

  * :ref:`Follow this link <tmslist_fnamecoil>`.


* **target**: *(3x1) list/array (Python/MATLAB) of floats*

  * **Desctiption**: Position of optimization target. Please see :ref:`coords_doc` for more information and how to obtain coordinates in SimNIBS
  * **Example**: *Python/MATLAB*

  .. code-block:: matlab

     opt.target = [-43.4, -20.7, 83.4];

  \ 

* **target_direction**: *None/[] or (3x1) list/array (Python/MATLAB) of floats, optional*

  * **Desctiption**: Direction of the electric field to be optimized. If :code:`None` (Python) or :code:`[]` (MATLAB), will optimize the elctric field magnitude. Default: optimize field magnitude
  * **Example**: *Python/MATLAB*

  .. code-block:: matlab

     opt.target_direction = [1, 0, 0]

  \ 


* **cond**: *list/array of COND structures (Python/MATLAB), optional*
   
  * :ref:`Follow this link <cond_attribute_doc>`.


* **anisotropy_type**: *'scalar', 'vn', 'dir' or 'mc', optional*

  * :ref:`Follow this link <anisotropy_type_attribute_doc>`.

* **aniso_maxratio**: *float, optional*

  * :ref:`Follow this link <aniso_maxratio_doc>`.

* **aniso_maxcond**: *float, optional*

  * :ref:`Follow this link <aniso_maxcond_doc>`.


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

  * **Description**: Reference position for the coil *y*-axis, with respect to the *pos* variable 
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

* **solver_options**: *str, optional*

  * **Description**: Options for the FEM matrix solver. Leave empty to use the CG+AMG solver, use a string to be passed to `PETSc <https://www.mcs.anl.gov/petsc/index.html>`_ or set :code:`solver_options='pardiso'` to use the MKL Pardiso direct solver
  * **Default**: Use the CG solver with the AMG preconditioner

* **method**: *'direct' or 'ADM', optional*

  * **Description**: Method to be used for the solution. Can be either a :code:`direct` method, which will perform electric field simulations at each coil position or the Auxiliary Dipole Method (ADM, `Gomez 2020 <https://doi.org/10.1101/2020.05.27.120022>`_), which uses reciprocity and the fast multipole method (FMM) to massively accelerate the optimization. The ADM method is only compatible with :file:`.ccd` coil files or :file:`.tcd` coil files that only contain dipole elements.
  * **Default**: :code:`direct`
