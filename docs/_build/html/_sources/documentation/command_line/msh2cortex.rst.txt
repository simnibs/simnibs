.. _msh2cortex_docs:

msh2cortex
===========

Description
-------------

Interpolates fields to a surface in the middle of the cortex, and optinally transforms the fields to the FsAverage template.

.. attention:: Not available for :ref:`headreco_docs` models ran with the :code:`--no-cat` option

Usage example
--------------

1. Open a simulation results folder, for example :file:`ernie/simu/` if you ran :file:`ernie_simu.mat`

2. Run

.. code-block:: bash

  msh2cortex -i ernie_TDCS_1_scalar.msh -m ../m2m_ernie/ -o subject_overlays -f fsavg_overlays

\

3. Freesufer overlays for the subject files will be created in the :file:`subject_overlays/` folder and for the average template in the :file:`fsavg_overlays folder/`. 

Further notes
----------------

* Type :code:`msh2cortex -h` for more information and options
* This tool equivalent to selecting the **interpolate to cortical surface** or the  **tranform to fsaverage space** options in the :ref:`GUI <sim_opt>` or selecting the **map_to_surf** or the **map_to_fsavg** to *true* in the :ref:`SESSION structure <session_doc>`.
