.. _tms_optimize:


TMS Optimization
==================

SimNIBS can help finding the best TMS position for stimulating a certain target. This is
done by searching coil positions in a grid around the target position and turning the
coil at various angles.

Basic Setting
--------------
Setting up a TMS optimization is similar to set-up a TMS simulation. In the most basic
settings, you need to select the head model, the coil being used and the target position.
The target position has to be given in :ref:`SimNIBS coordinates <coords_doc>` and can be
determined using the *nifti* volumes produced by :ref:`headreco_docs`, :ref:`mri2mesh_docs` or by using the :ref:`mni2subject_coords <mni2subject_coords_docs>` command line tool.

For accelerating the simulations, SimNIBS can use the MKL Pardiso direct solver. However, this
solver uses approximately three times more memory than the standard solver.


Python
''''''

.. literalinclude:: ../../../simnibs/examples/optimization/tms_optimization.py
   :language: python


MATLAB
''''''

.. literalinclude:: ../../../simnibs/examples/optimization/tms_optimization.m
   :language: matlab


This will first generate a grid of positions and start simulating. After it is done
simulating, SimNIBS will return with the position that causes the largest electric field
norm at the target.

The optimization will create the output file :file:`ernie_TMS_optimize_Magstim_70mm_Fig8_nii.msh`: Gmsh `.msh` file with the  optimized field and coil position


Refining the Search
--------------------
If you already have a good idea of where the coil should be located or oriented, you can
refine the search by precisely specifying the search region, search angles and resolution.

Python
''''''

.. literalinclude:: ../../../simnibs/examples/optimization/tms_opt_refined.py
   :language: python


MATLAB
''''''

.. literalinclude:: ../../../simnibs/examples/optimization/tms_opt_refined.m
   :language: matlab


Auxiliary Dipole Method (ADM)
---------------------------------------

To use the Auxiliary Dipole Method (ADM), simply use a :file:`.ccd` coil file and set the :code:`method = 'ADM'`

.. code-block:: python

  tms_opt.fnamecoil = 'Magstim_70mm_Fig8.ccd'
  tms_opt.method = 'ADM'

\



Acknowledgements
------------------

We would like to thank Ole Numssen and Konstantin Weise for the help in developing this
functionality, and Luis Gomez for contributing the code for the ADM optimization

Further Reading
------------------
Please see :ref:`tmsoptimize_doc` for a detail description of all TMS optimization options

References
------------

`Weise, K., Numssen, O., Thielscher, A., Hartwigsen, G., Knösche, T.R. (in review) A novel approach to localize cortical TMS effects. bioRxiv, 595603. <https://doi.org/10.1101/595603>`_

`Gomez, L. J., Dannhauer, M., & Peterchev, A. V. (2020). Fast computational optimization of TMS coil placement for individualized electric field targeting. bioRxiv <https://doi.org/10.1101/2020.05.27.120022>`_



