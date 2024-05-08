Python Library Documentation
=============================

SimNIBS main programming language is Python. Here you can see documentation for the main
modules used in SimNIBS.

.. currentmodule:: simnibs

I/O Functions
--------------

.. autosummary::
   :toctree: auto
   :nosignatures:

    read_msh
    write_msh
    read_freesurfer_surface
    write_freesurfer_surface
    read_gifti_surface
    read_curv
    write_curv
    read_stl
    write_geo_spheres
    write_geo_text



Mesh Classes
-------------

.. autosummary::
   :toctree: auto
   :nosignatures:

   Msh
   Nodes
   Elements
   NodeData
   ElementData

Transformations
----------------

.. autosummary::
   :toctree: auto
   :nosignatures:

   subject2mni_coords
   mni2subject_coords
   subject_atlas

Utilities
----------
.. autosummary::
   :toctree: auto
   :nosignatures:

    templates
    get_atlas
    SubjectFiles
	utils.region_of_interest.RegionOfInterest
    
Simulations
------------
.. autosummary::
   :toctree: auto
   :nosignatures:

    simulation.calc_B
    simulation.calc_fields
    simulation.tdcs
    simulation.tdcs_neumann
    simulation.tdcs_leadfield
    simulation.tms_dadt
    simulation.tms_coil
    simulation.tms_many_simulations
    simulation.electric_dipole