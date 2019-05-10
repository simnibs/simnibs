.. simnibs documentation master file, created by
   sphinx-quickstart on Wed Aug 22 14:44:31 2018.

SimNIBS 
=======

**DOCUMENTATION FOR THE 3.0 BETA VERSION**

**SimNIBS** is a free and open source software package for the Simulation of Non-invasive Brain Stimulation. It allows for realistic calculations of the electric field induced by transcranial magnetic stimulation (TMS) and transcranial electric stimulation (TES).


A SimNIBS workflow consists of three main parts:
   * :ref:`Generate high-quality head models <head_modeling_tutorial>` from MR images
   * Set up and run simulation in the :ref:`graphical user interface <gui_tutorial>`, :ref:`MATLAB or Python<scripting_tutorial>`
   * :ref:`Visualize simulation results <visualization_tutorial>`, and transform them to *FsAverage* or *MNI* space.

.. centered:: :ref:`Install SimNIBS <simnibs_installer>`

.. centered:: :ref:`Tutorial <tutorial>`

.. image:: images/simnibs_frontpage.png


SimNIBS is copyrighted |copy| by its :ref:`authors <contributors>` and licensed under :download:`GPL v3 <../LICENSE.txt>`.

.. |copy|   unicode:: U+000A9 .. COPYRIGHT SIGN

.. note:: When you publish results based on SimNIBS, please :ref:`cite our papers <publications>`
.. warning:: SimNIBS is a research tool. Clinical usage is not supported or advised. In particular, SimNIBS was not tested to give accurate results in the presence of pathological condition

News
======

  * **Version 3.0 beta is now available**. It is a major update of SimNIBS 
    
.. toctree::
   :maxdepth: 2
   :hidden:

   installation/installation
   tutorial/tutorial
   documentation/documentation
   contributors
   publications
   changelog
   internal/internal.rst

