Python Library Documentation
=============================

SimNIBS main programming language is Python. Here you can see documentation for the main
modules used in SimNIBS.

.. currentmodule:: simnibs.mesh_tools

I/O Functions
--------------

.. autosummary::
   :toctree: auto
   :nosignatures:

    mesh_io.read_msh
    mesh_io.write_msh
    mesh_io.read_freesurfer_surface
    mesh_io.write_freesurfer_surface
    mesh_io.read_gifti_surface
    mesh_io.read_curv
    mesh_io.write_curv
    mesh_io.read_stl
    mesh_io.write_geo_spheres
    mesh_io.write_geo_text



Mesh Classes
-------------

.. autosummary::
   :toctree: auto
   :nosignatures:

   mesh_io.Msh
   mesh_io.Nodes
   mesh_io.Elements
   mesh_io.NodeData
   mesh_io.ElementData

Transformations
----------------
.. currentmodule:: simnibs.utils.transformations

.. autosummary::
   :toctree: auto
   :nosignatures:

   subject2mni_coords
   mni2subject_coords
   subject_atlas

Utilities
----------
.. currentmodule:: simnibs.utils

.. autosummary::
   :toctree: auto
   :nosignatures:

    file_finder.templates
    file_finder.get_atlas
    file_finder.SubjectFiles
    
Simulations
------------
.. currentmodule:: simnibs.simulation

.. autosummary::
   :toctree: auto
   :nosignatures:

    biot_savart.calc_B
    fem.calc_fields
    fem.tdcs
    fem.tdcs_neumann
    fem.tdcs_leadfield
    fem.tms_dadt
    fem.tms_coil
    fem.tms_many_simulations
    fem.electric_dipole