.. simnibs documentation master file, created by
   sphinx-quickstart on Wed Aug 22 14:44:31 2018.

==========
 SimNIBS 4
==========


**SimNIBS 4** is a free and open source software package for the Simulation of Non-invasive Brain Stimulation. It allows for realistic calculations of the electric field induced by transcranial magnetic stimulation (TMS) and transcranial electric stimulation (TES).

A SimNIBS workflow consists of three main parts:

   * :ref:`Generate high-quality head models <head_modeling_tutorial>` from MR images.

   * Set up and run simulation using our fast FEM solvers in the :ref:`graphical user interface <gui_tutorial>`, or via :ref:`MATLAB or Python scripts <scripting_tutorial>`.

   * :ref:`Visualize simulation results <visualization_tutorial>`, and easily transform them to *FsAverage* or *MNI* space via in-built functions.

SimNIBS offers several advanced options such as :ref:`Optimizations of TMS coil positions <tms_optimize>`, or TES optimizations for :ref:`single and multiple targets <tdcs_opt>` and :ref:`brain network targeting <tdcs_distributed_opt>`.

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
    animation-duration: 18s;
    animation-delay: 0s;
    }
  .slideshow > div:nth-of-type(2) {
    animation-delay: 9s;
    }
  @keyframes fading-slideshow {
    0% { opacity: 0; }
    10% { opacity: 1; }
    50% { opacity: 1; }
    60% { opacity: 0; }
    }
  }
  </style>
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
		, 3000); /* time point for changing bgimg2 */
		setTimeout( function myRefresh2(){
		i=i+1;
		if (i >= shuffled.length) { i=0; }
		document.getElementById("bgimg1").style="background-image:url(_static/gallery/"+shuffled[i]+")"}
		, 12000);} /* time point for changing bgimg1 */
	, 18000);
  </script>
  </embed>

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

======
 News
======

  * **Version 4.0.0** is a major update that introduces our new head modeling approach :ref:`charm <head_modeling_tutorial>`. *Charm* replaces the previous methods *mri2mesh* and *headreco* which are discontinued. In contrast to the previous methods, *charm* does not have any external dependencies (matlab, freesurfer) and is ready-to-run directly after installation. Head models created by *charm* are more accurate for non-brain tissues and include additional tissue types such as large blood vessels and spongy bone.

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
