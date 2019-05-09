.. simnibs documentation master file, created by
   sphinx-quickstart on Wed Aug 22 14:44:31 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

SimNIBS: Simulation of Non Invasive Brain Stimulation
=========================================================

**DOCUMENTATION FOR THE 3.0 BETA VERSION**

**SimNIBS** is a free and open source software package for the Simulation of Non-invasive Brain Stimulation. It allows for realistic calculations of the electric field induced by transcranial magnetic stimulation (TMS) and transcranial electric stimulation (TES).


A SimNIBS workflow consists of three main parts:
   * Use :ref:`command line scripts <head_modeling_tutorial>` to generate high-quality meshes of the head from MR images, or use our **example data set**.
   * Use the :ref:`graphical user interface <gui_tutorial>` or :ref:`scripts <scripting_tutorial>` to set up and run simulations.
   * :ref:`Visualize simulation results <visualization_tutorial>` in *gmsh* or *MATLAB*, transform them to *FsAverage* or *MNI* space, or load them in *MATLAB* or *Python*.

.. image:: images/simnibs_frontpage.png

:ref:`INSTALL SIMNIBS <installation>`


:ref:`TUTORIAL <tutorial>`

SimNIBS is copyrighted |copy| by its :ref:`authors <contributors>` and licensed under :download:`GPL v3 <../LICENSE.txt>`.

.. |copy|   unicode:: U+000A9 .. COPYRIGHT SIGN


.. warning:: SimNIBS is a research tool. Clinical usage is not supported or advised. In particular, SimNIBS was not tested to give accurate results in the presence of pathological condition

.. note:: When you publish results based on SimNIBS 2.0 or 2.1, or parts of it, please cite: `Thielscher, A., Antunes, A. and Saturnino, G.B. (2015), Field modeling for transcranial magnetic stimulation: a useful tool to understand the physiological effects of TMS? IEEE EMBS 2015, Milano, Italy <http://dx.doi.org/10.1109/EMBC.2015.7318340>`_

News
======

  * **Version 3.0 beta is now available**. It is a major update of SimNIBS 
    
.. toctree::
   :maxdepth: 2
   :hidden:

   installation
   tutorial/tutorial
   documentation/documentation
   contributors
   publications
   changelog
   internal/internal.rst

