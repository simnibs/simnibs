.. _group_tutorial:

Group Analysis
===============


SimNIBS supports group analysis of electric fields in both *FsAverage* and *MNI* spaces.

In this example, we will simulate the average electric field for a simple tDCS montage across five subjects using the *FsAverage* space.

Example Dataset
----------------


`Click here to download the example dataset for group analyses <https://osf.io/egshq/download>`_.

This example dataset is composed of a subgroup of a cohort available at `OpenfMRI <https://openneuro.org/datasets/ds000171>`_. The data was processed in SimNIBS 2.1 using :ref:`headreco_docs`. For more information, please see the `OSF repository <https://osf.io/ah5eu/>`_ and `Saturnino et al., 2018 <https://doi.org/10.1101/500314>`_.


In this example, we have five subjects, named :file:`sub01`, :file:`sub09`,
:file:`sub10`, :file:`sub12`, and :file:`sub15`.


Objective
---------

Suppose we apply tDCS on the five subjects using an anode over C3 and a cathode over AF4.
We want to obtain the mean and standard deviation of the normal component of the electric field (that is, the incoming or outgoing electric field component) in the middle gray matter surface across all subjects.

.. figure:: ../../images/tutorial_normal.png

   Norm (top) and Norma (bottom) of the electric field, from `Antonenko et al. 2019 <https://doi.org/10.1016/j.brs.2019.03.072>`_

\


Set up and Run Simulations
---------------------------

There are several ways to set-up and run Simulations in SimNIBS


GUI
''''
Set-up the simulation for each subject in the :ref:`Graphical User Interface <gui_tutorial>`. In this case, remember to tick the *Transform to FsAverage* option in the :ref:`sim_opt`. (under *Edit* -> *Simulation Options*)

Python
'''''''
Write a *Python* script. In this case, remember to set *map_to_fsavg* to *True* in the :ref:`session_doc` structure. See :ref:`scripting_tutorial` for more information.

.. literalinclude:: ../../../simnibs/examples/simulations/run_simulations_group.py
   :language: python

MATLAB
''''''''
Write a *MATLAB* script. In this case, remember to set *map_to_fsavg* to *True* in the :ref:`session_doc` structure. See :ref:`scripting_tutorial` for more information.


.. literalinclude:: ../../../simnibs/examples/simulations/run_simulations_group.m
   :language: matlab




Calculate Mean 
----------------

When the simulations are over, we need to collect their results to calculate an average. In SimNIBS, we can do it either in Python or MATLAB.
Please notice that, while for setting up simulations Python and MATLAB share a similar interface, in post-processing the interfaces can be very different.

Python
''''''''

.. literalinclude:: ../../../simnibs/examples/analysis/group_average_fsavg.py
   :language: python


MATLAB
'''''''
.. literalinclude:: ../../../simnibs/examples/analysis/group_average_fsavg.m
   :language: matlab




Further Reading
---------------
For more information on group analysis, please see our `SimNIBS 2.1 tutorial paper <https://doi.org/10.1101/500314>`_ and `Bungert et al. 2017 <https://doi.org/10.1093/cercor/bhw292A>`_.

