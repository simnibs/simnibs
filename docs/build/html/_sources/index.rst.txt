.. simnibs documentation master file, created by
   sphinx-quickstart on Wed Aug 22 14:44:31 2018.

SimNIBS 
=======


**SimNIBS** is a free and open source software package for the Simulation of Non-invasive Brain Stimulation. It allows for realistic calculations of the electric field induced by transcranial magnetic stimulation (TMS) and transcranial electric stimulation (TES).

A SimNIBS workflow consists of three main parts:
   * :ref:`Generate high-quality head models <head_modeling_tutorial>` from MR images
   * Set up and run simulation in the :ref:`graphical user interface <gui_tutorial>`, :ref:`MATLAB or Python<scripting_tutorial>`
   * :ref:`Visualize simulation results <visualization_tutorial>`, and transform them to *FsAverage* or *MNI* space.

|

.. image:: images/simnibs_frontpage.png

|

.. raw:: html

  <style>
  * {
    box-sizing: border-box;
  }
  div {
    text-align: left;
  }
  .column {
    float: left;
    width: 50%;
    padding: 10px;
  }
  
  .row:after {
    content: "";
    display: table;
    clear: both;
  }
  </style>
  <div class="row">
  <div class="column">
  SimNIBS is being developed at the <a href="https://www.drcmr.dk">Danish Research Centre for Magnetic Resonance (DRCMR)</a> and the <a href="https://www.dtu.dk/english">Technical University of Denmark (DTU)</a>, in collaboration with <a href="contributors.html" >external partners</a>
  </div>
  <div class="column">
  <a href="https://www.drcmr.dk"><img src="_static/DRCMR_logo.png" alt="DRCMR", height="70"></a>
  <a href="https://www.dtu.dk/english"><img src="_static/DTU_logo.png" alt="DTU", height="70"></a>
  </div>
  </div>


SimNIBS is copyrighted |copy| by its :ref:`authors <contributors>` and licensed under :download:`GPL v3 <../LICENSE.txt>`.


.. note:: When you publish results based on SimNIBS, please cite `Thielscher, A., Antunes, A. and Saturnino, G.B. (2015), Field modeling for transcranial magnetic stimulation: a useful tool to understand the physiological effects of TMS? IEEE EMBS 2015, Milano, Italy <http://dx.doi.org/10.1109/EMBC.2015.7318340>`_ 

  Please see :ref:`our publications <publications>` for details about the SimNIBS modules.

.. warning:: SimNIBS is a research tool. Clinical usage is not supported or advised. In particular, SimNIBS was not tested to give accurate results in the presence of pathological condition

.. |copy| unicode:: U+000A9

News
======

  * **Version 3.2 is now avaliable**  featuring :ref:`TMS optimization  using the ADM method <tms_optimize>` (contributed by Luis Gomez), TES :ref:`field strength <tdcs_opt>` and :ref:`brain network <tdcs_distributed_opt>` optimization and magnetic field calculations for MRCDI and MREIT 
  * **Version 3.1 is now avaliable**. This version features :ref:`TDCS electrode optimization <tdcs_opt>`, :ref:`TMS coil position optimization <tms_optimize>` (contributed by Konstantin Weise and Ole Numssen) and even faster FEM solvers for leadfield calculations
  * **Version 3.0 is now available**. It is a major update of SimNIBS. Expect much faster simulations and easier visualizations!
    
.. toctree::
   :maxdepth: 2
   :hidden:

   installation/installation
   dataset
   tutorial/tutorial
   documentation/documentation
   faq
   contributors
   publications
   changelog
   contact
