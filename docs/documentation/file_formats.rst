.. _file_formats:

File Formats in SimNIBS
=========================

Head Models
------------

SimNIBS' head models are **meshes** (*.msh* files).
This means that the head head is represented as a set of *Nodes* and *Elements*.

    * The *Nodes* are points located in the 3-dimensional volume.
    * The *Elements* are triangles and tetrahedra. They are defined using 3 nodes (triangles) or 4 nodes (tetrahedra).

This type of format is highly advantageous for Finite Element (FEM) calculations, especially for complex geometries such as the human head.

Head Meshes are stored in **binary gmsh version 2 format**, as described in the `Gmsh documentation <http://gmsh.info/doc/texinfo/gmsh.html#MSH-file-format-version-2>`_.
SimNIBS can read arbitrary mesh files.
But it can't write files with elements other then first order triangles and tetrahedra.

SimNIBS offers the *Python* :func:`simnibs.msh.read_msh` function and the
*MATLAB* function *mesh_load_gmsh4* to read *.msh* files.

TMS Coils
---------
Definitions of TMS coils can be stored in a number of different file formats:

  * Dipole files (.ccd) which store the position and magnitude of dipoles 
  * Files in the NiftI format (.nii | .nii.gz) which store sampled values of the magnetic vector potential (A-field) on a grid
  * A combined TMS coil definition file format (.tcd) which can store dipole elements, line segments and a sampled magnetic vector potential (A-field) on a grid. This format also allows for flexible or moveable coil components.

TMS Coil Definition File Format (.tcd)
''''''''''''''''''''''''''''''''''''''
The TMS Coil Definition File Format (.tcd) is a JSON-based format that defines a single TMS coil. 
The coil file can be loaded into SimNIBS and used to calculate the magnetic vector potential (A-field) for a given rate of change
of the coil current and the magnetic field (B-field) for a given coil current. The following text gives an overview of the format.
Details are provided in the SimNIBS documentation as well as in the example scripts for custom coil creation. 

.. figure:: ../images/tcd_json.svg

    The TMS Coil Definition File Format. 
    Each sub list (green) can have an arbitrary number of items. 
    A) Each coil can have a general coil casing (or “coilModel”) not associated with any coil element. 
    B) Each coil element can have a casing. C) Each coil element must be connected to one stimulator. 
    D) Each coil element can have an arbitrary number of deformations including no deformations. 
    E) Each deformation uses a deformation range that defines the valid values of the deformation. 
    F) Each coil can have an arbitrary number of self-intersection groups.

The file format holds basic information about the coil, the name, the brand, and the version of the coil file.
Its main part is a list of stimulating coil elements (“coilElementList”) that is linked to a list of deformations ("deformList") and a list
of geometric models of the sub-casing (“coilModels”). 
These three components allow for definitions of static coils and coils with one or multiple movable stimulating elements.
The following list gives a brief overview of the main components:

    * *coilElementList*

      Each stimulating coil element can use one of three types of stimulating elements:
      
      * Magnetic dipoles, which are described by their location (in mm), direction, and magnitude (in Am²).
      * Line current elements, which are defined by their starting position (in mm), direction, and magnitude (in mm). 
      * Magnetic vector potential field (A-field), sampled on a grid, described by the A-field value at each grid location (in Tm) and an affine transformation matrix (in mm).
      
      Each stimulating coil element has a name, is connected to one stimulator (from “stimulatorList”)
      (arrow “C” in the image), has a casing (from “coilModels”; Image arrow B),
      and is associated with a list of deformations (from “deformList”) 
      that are applied to this coil element (Image arrow D).
      The points, values and data lists are stored as plain text or as Base64 encoded and compressed binary.
      The deformations are applied in the order that they are stored in.
      Multiple stimulating coil elements can share the same deformation,
      enabling the common rotation of several coil elements around shared axes.
    * *stimulatorList*

      A list of TMS stimulators is stored to connect different stimulating coil elements to the same stimulator, which can be used in multi-stimulator coil settings. The TMS stimulators are described by their name, brand, maximum dI/dt (in A/s), and a list of waveforms.
    * *deformList*

      A deformation is either defined as a translation in the x, y, or z direction (in mm) or as a rotation around an axis (in degrees) defined by two points (in mm). Its value is defined and bound by a deformation range (Image arrow E).
    * *deformRangeList*

      The limits for deformations are defined as deformation ranges. Multiple deformations, for example, rotations around different axes, can share the same deformation range. These ranges are defined by a minimum and maximum value, and an initial value between the minimum and the maximum.
    * *coilModels*

      For visualization and optimization purposes, the file format contains a list of coil casing triangulated surfaces. Each coil casing surface has a point list (in mm) and a list of faces which are indices into the point list. Additionally, a list of points ("minDistancePoints" in mm) representing the parts of the coil that are supposed to be close to the head can be stored. These are then used for evaluation of the coil-skin distance costs during optimization. The point and faces lists are stored as plain text or as Base64 encoded and compressed binary.
    * *selfIntersectionTest*  

      A list of groups of stimulating coil element indices which will be tested for self-intersections in cases of optimization. 

The file format also gives the possibility to store a global coil casing ("coilCasing") that is not associated with a specific stimulating coil element (Image arrow A). In addition, to speed up repeated FEM evaluations for varying coil positions and deformations, the A-fields of the stimulating coil elements can be sampled on 3D grids during preparation. They can then be computationally efficient evaluated for varying coil configurations using simple linear transformations. The file format supports this by giving the option to store information about the axis limits ("limits" in mm) and resolution ("resolution" in mm) of the 3D grids.

Examples on how to create TMS coil definitions in the tcd file format can be found under :file:`simnibs/examples/coils`.
The JSON schema to validate tcd files can be found under :file:`simnibs/resources/coil_models/coil_model.schema.json`.

Simulation Results
----------------------

By default, Simulation results are also stored in *Gmsh* format.
There are 2 types of fields used in SimNIBS:

    * *NodeData* is defined at each nodes. By default, only the electric potential "v" is stored as *NodeData*

    * *ElementData* is defined for each element. This is the format of choice for the electric field, current density and their respective magnitudes.

The choice of format for each field is due to how the Finite Element Method works.
After the FEM calculations, we obtain values of "v" at each node.
To obtain the electric field and the current density, we take the gradient of the potential.
However, the gradient operation is defined element-wise, and not node-wise.

Together with the *.msh* files, we often also save *.opt* files to facilitate visualization.

Surface Fields
---------------

When fields are mapped to the middle gray matter surface, either on the subject or on the  *FsAverage*, it saves results as a *FreeSurfer* *.curv* file, which contains values for each point in the surface. SimNIBS has the *mesh_load_fsresults* *MATLAB* function and the :func:`simnibs.msh.read_curv` *Python* function to load this kind of file.

The surfaces themselves are stored as GIFTI files in :file:`m2m_{subID}/surfaces/`. They can be read using *mesh_load_fssurf* in *MATLAB* and :func:`simnibs.msh.read_gifti_surface` in *Python*.


Volumes
--------

Fields mapped to subject or MNI volumes are stored in NiftI format.

HDF5
----

SimNIBS uses `HDF5 <https://www.hdfgroup.org/solutions/hdf5/>`_ to store large data sets, such as for uncertainty quantification (UQ) and leadfields.
The HDF5 format is hierarchical, meaning that is acts almost as a folder structure.

