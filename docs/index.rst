.. simnibs documentation master file, created by
   sphinx-quickstart on Wed Aug 22 14:44:31 2018.

=============
 SimNIBS 4.5
=============


**SimNIBS v4.5** is an open source software package for the Simulation of Non-invasive Brain Stimulation. It allows for realistic calculations of the electric field induced by transcranial magnetic stimulation (TMS) and transcranial electric stimulation (TES).

A basic SimNIBS workflow consists of three main parts:

   * :ref:`Generate high-quality head models <head_modeling_tutorial>` from MR images.

   * Set up and run simulation using our fast FEM solvers in the :ref:`graphical user interface <gui_tutorial>`, or via :ref:`MATLAB or Python scripts <scripting_tutorial>`.

   * :ref:`Visualize simulation results <visualization_tutorial>`, and easily transform them to *FsAverage* or *MNI* space via in-built functions.

SimNIBS offers several advanced options such as :ref:`Optimizations of TMS coil positions <overview_tms_opt>` or several :ref:`TES optimization methods <overview_tes_opt>`. Also a number of :ref:`external open-source software packages and workflows <external_packages>` build upon SimNIBS.

.. raw:: html

  <embed>
  <style>
  .slideshow {
    position: relative;
    width: 100%;
    height: 30vh;
    background-color: #fff;
    }
  .slideshow > div {
    width: 100%;
    height: 100%;
    background-size: contain;
    background-repeat: no-repeat;
    background-position: center;
    background-image: url();
    position: absolute;
    opacity: 0;
    animation-name: fading-slideshow;
    animation-iteration-count: infinite;
    animation-duration: 12s;
    animation-delay: 0s;
    }
  .slideshow > div:nth-of-type(2) {
    animation-delay: 6s;
    }
  @keyframes fading-slideshow {
    0% { opacity: 0; }
    10% { opacity: 1; }
    50% { opacity: 1; }
    60% { opacity: 0; }
    }
  }
  </style>
  <a href="gallery.html">
  <div class="slideshow">
    <div id="bgimg1";></div>
    <div id="bgimg2";></div>
  </div>
  <script type="text/javascript" src=_static/gallery/list.js></script>
  <script type="text/javascript">
    let shuffled = filelist
		.map(value => ({ value, sort: Math.random() }))
		.sort((a, b) => a.sort - b.sort)
		.map(({ value }) => value)	
	document.getElementById("bgimg1").style="background-image:url(_static/gallery/"+filelist[0]+")";  
	document.getElementById("bgimg2").style="background-image:url(_static/gallery/"+shuffled[0]+")"
	i=0
	setTimeout( function myRefresh2(){
		i=i+1;
		if (i >= shuffled.length) { i=0; }
		document.getElementById("bgimg1").style="background-image:url(_static/gallery/"+shuffled[i]+")"}
		, 12000); /* time point for changing bgimg1 */
	setInterval( function myRefresh(){
		setTimeout( function myRefresh1(){
			i=i+1;
			if (i >= shuffled.length) { i=0; }
			document.getElementById("bgimg2").style="background-image:url(_static/gallery/"+shuffled[i]+")"}
		, 6000); /* time point for changing bgimg2 */
		setTimeout( function myRefresh2(){
		i=i+1;
		if (i >= shuffled.length) { i=0; }
		document.getElementById("bgimg1").style="background-image:url(_static/gallery/"+shuffled[i]+")"}
		, 12000);} /* time point for changing bgimg1 */
	, 12000);
  </script>
  </embed>
  </a href="gallery.html">

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

It builds upon 3rd party open-source and free code and binaries from :ref:`several projects and research groups. <3rd_party>`

.. note:: When you publish results based on SimNIBS, please cite `Thielscher, A., Antunes, A. and Saturnino, G.B. (2015), Field modeling for transcranial magnetic stimulation: a useful tool to understand the physiological effects of TMS? IEEE EMBS 2015, Milano, Italy <http://dx.doi.org/10.1109/EMBC.2015.7318340>`_

  Please see :ref:`our publications <publications>` for details about the SimNIBS modules.

.. warning:: SimNIBS is a research tool. Clinical usage is not supported or advised. In particular, SimNIBS was not tested to give accurate results in the presence of pathological condition

.. |copy| unicode:: U+000A9

======
 News
======

**Version 4.5.0** includes:

* Exciting new optimization methods for :ref:`TMS <tms_flex_opt>`  and :ref:`TES <tes_flex_opt>`.
	* :ref:`Optimization of TMS coil positions also for bent and flexible coils <tms_flex_opt>`, thereby systematically avoiding intersections of the coil with the head.
	* :ref:`Leadfield-free optimization of TES montages <tes_flex_opt>`, including those with rectangular electrodes, center-surround montages, temporal interference stimulation and electrode arrays for tumor treating field therapies.
* Several new TMS coil models and new :ref:`datasets <dataset>`.
* New format for TMS coil models (.tcd) that supports flexible and multi-element coils and simplifies the creation of custom coil models (example scripts are provided).
* Tutorial for calculating :ref:`EEG leadfields <eeg_leadfields>` with SimNIBS for use in `FieldTrip <https://www.fieldtriptoolbox.org/>`_ and `MNE-Python <https://mne.tools/stable/index.html>`_ .
* JupyterLab to make SimNIBS scripting in Python more straightforward.

Please see details in the :ref:`changelog <changelog>`.

.. toctree::
   :maxdepth: 2
   :hidden:

   installation/installation
   dataset
   tutorial/tutorial
   documentation/documentation
   faq
   gallery
   contributors
   publications
   changelog
   contact
   external_packages
