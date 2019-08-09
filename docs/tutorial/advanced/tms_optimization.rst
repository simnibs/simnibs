.. _tms_optimize:


TMS Optimization
==================

SimNIBS can help finding the best TMS position for stimulating a certain target. This is
done by searching coil positions in a grid around the target position, at various angles.

Basic Setting
--------------
Setting up a TMS optimization is similar to set-up a TMS simulation. In the most basic
settings, you need to select the head model, the coil being used and the target position.
The target position are needed as :ref:`SimNIBS coordinates <coords_doc>` and can be
determined using the *nifti* volumes produced by :ref:`headreco_docs`, :ref:`mri2mesh_docs` or by using the :ref:`mni2subject_coords <mni2subject_coords_docs>` command line tool.

Python
''''''

.. literalinclude:: ../../../simnibs/examples/tms_optimization.py
   :language: python


MATLAB
''''''

.. literalinclude:: ../../../simnibs/examples/tms_optimization.m
   :language: matlab


This will first generate a grid of positions and start simulating. After it is done
simulating, SimNIBS will return with the position that causes the largest electric field
norm at the target.

For accelerating the repeated simulations, SimNIBS will use the MKL Pardiso direct solver. This
means that TMS optimization has therefore a larger memory footprint than regular simulations.


Refining the Search
--------------------
If you already have a good idea of where the coil should be located or oriented, you can
refine the search by precisely specifying the search region and resolution.

Python
''''''

.. literalinclude:: ../../../simnibs/examples/tms_optimization_refined.py
   :language: python


MATLAB
''''''

.. literalinclude:: ../../../simnibs/examples/tms_optimization_refined.m
   :language: matlab


Acknolowedgements
------------------

We would like to thank Ole Numssen and Konstantin Weise for the help in developing this
functionality
