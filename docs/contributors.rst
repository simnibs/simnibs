.. _contributors:

Contributors and Funding
=========================

Current Contributors
---------------------
* Axel Thielscher
   * Project coordination
   * Debugging and maintenance across the complete code.
   * Meshing part of :ref:`charm_docs` (4.0)
   * :ref:`dwi2cond_docs` (2.1)
   * Example dataset *ernie* (2.1)

* Oula Puonti
   * Main author of :ref:`charm_docs` and of simnibs 4.0 in general (4.0)
   * Debugging and maintenance across the complete code.
   * New MNI head mesh (2.1)

* Kristoffer H. Madsen
   * Debugging and maintenance across the complete code.
   * Cross-platform building (4.0)
   * new html-based viewer (4.0)
   * ccd to nifti conversion (4.0)
   * MNI transformation (2.1)
   * Fast I/O for gmsh-meshes in python (2.0)
   * Pre-calculated A-fields for TMS in the pipeline using nifti volume files (2.0)

* Jesper D. Nielsen
   * Contributions to the segmentation functions in :ref:`charm_docs` (4.0)
   * Main author of the *headreco* pipeline, now discontinued (2.1)

* Fang Cao
   * Code testing and updating to python 3.9 (4.0)

* Konstantin Weise
   * First version of the TMS optimization (3.1)
   * UQ functionality (3.0)

* Thomas Knoesche
   * Help with the UQ functionality (3.0)

* Ole Numssen
   * First version of the TMS optimization (3.1)

Former Contributors
---------------------

* Guilherme B. Saturnino
   * Main contributor to many SimNIBS features: 
	   * TES optimiaztion algorithms (3.1, 3.2)
	   * Installation procedure (3.0, 3.2)
	   * Documentation (3.0, 3.1 and 3.2)
	   * Fast FEM code (3.0)
	   * (together with K. Weise) UQ functionality (3.0)
	   * New python core (2.1)
	   * GUI (2.0)
	   * Electrode modeling for TES (2.0)
	   * Bug-fixing
	   * Meshing part of :ref:`charm_docs` (4.0)

* Hassan Yazdanian and Kim Knudsen
   * Magnetic Field Calculations (3.2)

* Luis J. Gomez, Moritz Dannhauer, and Angel V. Peterchev; Duke University, Durham, North Carolina, U.S.A.
   * Auxiliary Dipole Method (ADM) TMS optimization (3.2)

* Andre Antunes
   * Main contributor to the FEM pipeline in SimNIBS 2.0
   * Implementation of a range of post-processing programs in SimNIBS 2.0
   
* Andreas Bungert
   * Testing of the new FEM pipeline in SimNIBS 2.0

* Alex Opitz
   * Implementation of first diffusion-to-conductivity mapping approach
   * Co-contributor to many other parts in SimNIBS 1.0
   * Testing and validation of the new FEM calculations for tDCS in SimNIBS 2.0

* Mirko Windhoff
   * Main contributor to SimNIBS 1.0
   
Acknowledgements 
-----------------
SimNIBS integrates free software for Neuroimaging, computer graphics
and FEM calculations into one coherent pipeline:

* :ref:`charm_docs` uses `Samseg (Oula Puonti, Koen Van Leemput) from FreeSurfer <https://surfer.nmr.mgh.harvard.edu/fswiki/Samseg>`_ as segmentation backend, `CGAL <https://www.cgal.org/>`_ for meshing, and also a modified version of `MeshFix <http://code.google.com/p/meshfix/>`_ by `Marco Attene <https://www.cnr.it/en/people/marco.attene>`_, functions from `CAT12 <http://dbm.neuro.uni-jena.de/cat/>`_ from the `Structural Brain Mapping Group (University of Jena) <http://www.neuro.uni-jena.de/>`_
* We use `Gmsh <http://geuz.org/gmsh/>`_ by `Christophe Geuzaine 
  <http://www.montefiore.ulg.ac.be/~geuzaine/>`_ and `Jean-François Remacle <http://perso.uclouvain.be/jean-francois.remacle/>`_ for visualization.
* :ref:`dwi2cond_docs` uses `FSL <http://www.fmrib.ox.ac.uk/fsl/>`_ from the `FMRIB Center (Oxford University) <http://www.fmrib.ox.ac.uk/>`_.
* FEM code introduced in version 3.0 relies on `PETSc <https://www.mcs.anl.gov/petsc/>`_ and `Hypre
  <https://computation.llnl.gov/projects/hypre-scalable-linear-solvers-multigrid-methods/software>`_.
* Thanks for **Konstantin Weise** for giving us access to an early version of his software, `pygpc <https://github.com/konstantinweise/pygpc>`_.

* Thanks for **Luis Gomez** for providing a python version of his `ADM TMS optimization package <https://github.com/luisgo/Auxiliary_dipole_method>`_.

Institutions
---------------

* Versions 2.1, 3 and 4 have been developed at the `Danish Research Center for Magnetic Resonance <http://www.drcmr.dk>`_ (Copenhagen, Denmark) and the `Technical University of Denmark <http://www.dtu.dk/english>`_ (Kgs Lyngby, Denmark), in collaboration with external partners.
* Version 1.0 was created at the `Max-Planck Institute for Biological Cybernetics <http://www.kyb.tuebingen.mpg.de>`_ (Tübingen, Germany).
* Version 2.0 was developed in all three institutions

Funding Sources
-----------------

We would like to thank our funding sources

.. centered::  |lundbeck|_ |novo|_ |sdc|_ |stiped|_ |if|_ |nimh|_


.. |lundbeck| image:: ./images/lundbeckfonden.png
   :height: 50
.. _lundbeck: https://www.lundbeckfonden.com/en/

.. |novo| image:: ./images/novonordiskfonden.png
   :height: 50
.. _novo: https://novonordiskfonden.dk/en/

.. |sdc| image:: ./images/sdc.png
   :height: 50
.. _sdc: http://sdc.university/

.. |stiped| image:: ./images/stiped.png
   :height: 50
.. _stiped: http://www.stiped.eu/home/

.. |if| image:: ./images/innovationsfonden.png
   :height: 50
.. _if: https://innovationsfonden.dk/en

.. |nimh| image:: ./images/NIH-NIMH-logo-new.png
   :height: 50
.. _nimh: https://www.nimh.nih.gov/index.shtml


