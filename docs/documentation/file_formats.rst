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

