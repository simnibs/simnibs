.. _contributors:

Contributors
===============

Current Contributors
---------------------
* Axel Thielscher
   * Project coordination
   * Debugging and maintenance of mri2mesh, headreco and meshfix
   * Implementation of new diffusion-to-conductivity mapping script (dwi2cond)
   * Creation of the new example dataset *ernie*

* Guilherme B. Saturnino
   * Implementation of the new FEM code in SimNIBS 3.0
   * Main contributor to the new python core of SimNIBS 2.1
   * Main contributor to the GUI and the electrode modeling for tDCS in 2.0
   * Implementation of the new installation procedure in 2.0
   * Bug-fixing

* Jesper D. Nielsen
   * Main author of the new **headreco** script in 2.1

* Kristoffer H. Madsen
   * Implementation of prototype version for MNI transformation in 2.1
   * Implementation of fast I/O for gmsh-meshes in python in SimNIBS 2.0
   * Integration of pre-calculated A-fields for TMS in the pipeline using nifti volume files in 2.0

* Oula Puonti
   * Help with **headreco**
   * Creation of the new MNI head mesh

Former Contributors
---------------------
* Mirko Windhoff
   * Main contributor to SimNIBS 1.0

* Alex Opitz
   * Implementation of first diffusion-to-conductivity mapping approach
   * Co-contributor to many other parts of the first version
   * Testing and validation of the new FEM calculations for tDCS in SimNIBS 2.0

* Andre Antunes
   * Main contributor to the FEM pipeline in SimNIBS 2.0
   * Implementation of a range of post-processing programs in SimNIBS 2.0

* Andreas Bungert
   * Testing of the new FEM pipeline in SimNIBS 2.0

Acknowledgements 
-----------------
SimNIBS integrates free software for Neuroimaging, computer graphics
and FEM calculations into one coherent pipeline:

* **mri2mesh** uses `FreeSurfer <http://surfer.nmr.mgh.harvard.edu/>`_ from the
  `Athinoula A. Martinos Center for Biomedical Imaging <http://www.nmr.mgh.harvard.edu/martinos/flashHome.php>`_ and `FSL
  <http://www.fmrib.ox.ac.uk/fsl/>`_ from the `FMRIB Center (Oxford University) <http://www.fmrib.ox.ac.uk/>`_.
* **headreco** uses `SPM12 <https://www.fil.ion.ucl.ac.uk/spm/software/spm12/>`_ from the
  Wellcome Trust Centre for Neuroimaging (UC London) and `CAT12 <http://dbm.neuro.uni-jena.de/cat/>`_ from the `Structural Brain Mapping Group (University of Jena) <http://www.neuro.uni-jena.de/>`_
* Both use a modified version of `MeshFix <http://code.google.com/p/meshfix/>`_ by `Marco Attene <http://pers.ge.imati.cnr.it/attene/PersonalPage/attene.html>`_. Many thanks go to Marco for releasing MeshFix as open source and for his support when extending it!
* We heavily use `Gmsh <http://geuz.org/gmsh/>`_ by `Christophe Geuzaine
  <http://www.montefiore.ulg.ac.be/~geuzaine/>`_ and `Jean-François Remacle <http://perso.uclouvain.be/jean-francois.remacle/>`_
* FEM code introduced in version 3.0 relies on `PETSc <https://www.mcs.anl.gov/petsc/>`_ and `Hypre
  <https://computation.llnl.gov/projects/hypre-scalable-linear-solvers-multigrid-methods/software>`_.
* Thanks for Konstantin Weise for giving us access to an early version of his software, `pygpc <https://github.com/konstantinweise/pygpc>`_.

Institutions
---------------

* Versions 2.1 and 3.0 have been developed at the `Danish Research Center for Magnetic Resonance <http://www.drcmr.dk>`_ (Copenhagen, Denmark) and the `Technical University of Denmark <http://www.dtu.dk/english>`_ (Kgs Lyngby, Denmark).
* Version 1.0 was created at the `Max-Planck Institute for Biological Cybernetics <http://www.kyb.tuebingen.mpg.de>`_ (Tübingen, Germany).
* Version 2.0 was developed in all three institutions


