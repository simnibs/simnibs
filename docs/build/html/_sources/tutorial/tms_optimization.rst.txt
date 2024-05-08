.. _tms_optimize:


TMS Optimization based on Grid Search or the Auxiliary Dipole Method (ADM)
==========================================================================

SimNIBS can help finding the best TMS position for stimulating a certain target. This is
done by searching coil positions in a grid around the target position and turning the
coil at various angles.

.. note::  

   When using this feature in a publication, please cite either

   `Weise, K., Numssen, O., Thielscher, A., Hartwigsen, G., & Knösche, T. R. (2020). A novel approach to localize cortical TMS effects. Neuroimage, 209, 116486. <https://doi.org/10.1016/j.neuroimage.2019.116486>`_
   
   or

   `Gomez, L. J., Dannhauer, M., & Peterchev, A. V. (2020). Fast computational optimization of TMS coil placement for individualized electric field targeting. NeuroImage 2021; 228: 117696. <https://doi.org/10.1016/j.neuroimage.2020.117696>`_
   
   in case you use the ADM method. The ADM code is distributed under the conditions below:


   Authors of auxiliary dipole method (ADM) and code for determining rapidly the optimum coil position and orientation: Luis J. Gomez, Moritz Dannhauer, and Angel V. Peterchev; Duke University, Durham, North Carolina, U.S.A.

   The development of the Duke ADM algorithm and code were supported by the National Institute of Mental Health and the National Institute on Aging of the National Institutes of Health under Award Numbers K99MH120046, RF1MH114268, RF1MH114253, and U01AG050618. The content is solely the responsibility of the authors and does not necessarily represent the official views of the National Institutes of Health.
   
   The copyrights of this software are owned by Duke University. As such, two licenses to this software are offered:

      * An open source license under the GNU General Public License (GPL-v2.0) (https://opensource.org/licenses/gpl-2.0.php).
      * A custom license with Duke University, for use without the GPL-v2.0 restrictions. 

   As a recipient of this software, you may choose which license to receive the code under. Outside contributions to the Duke owned code base cannot be accepted unless the contributor transfers the copyright to those changes over to Duke University.

   To enter a license agreement without the GPL-v2.0 restrictions, please contact the Digital Innovations department at Duke Office of Licensing and Ventures (https://olv.duke.edu/software/) at olvquestions@duke.edu with reference to “OLV File No. 7148” in your email. 

   Please note that this software is distributed AS IS, WITHOUT ANY WARRANTY; and without the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

 
\

Basic Setting
--------------
Setting up a TMS optimization is similar to setting-up a TMS simulation. In the most basic
setting, you need to select the head model, the coil being used and the target position.
The target position has to be given in :ref:`SimNIBS coordinates <coords_doc>` and can be
determined using the *nifti* volumes produced by :ref:`headreco_docs`, :ref:`mri2mesh_docs` or by using the :ref:`mni2subject_coords <mni2subject_coords_docs>` command line tool.

For accelerating the simulations, SimNIBS can use the MKL Pardiso direct solver. However, this
solver uses approximately three times more memory than the standard solver.


Python
''''''

.. literalinclude:: ../../simnibs/examples/optimization/tms_optimization.py
   :language: python


MATLAB
''''''

.. literalinclude:: ../../simnibs/examples/optimization/tms_optimization.m
   :language: matlab


This will first generate a grid of positions and start simulating. After it is done
simulating, SimNIBS will return with the position that causes the largest electric field
magnitude at the target.

The optimization will create the Gmsh output file :file:`ernie_TMS_optimize_Magstim_70mm_Fig8_nii.msh` with the optimized field and coil position


Refining the Search
--------------------
If you already have a good idea of where and how the coil should be located or oriented, you can
refine the search by precisely specifying the search region, search angles and resolution.

Python
''''''

.. literalinclude:: ../../simnibs/examples/optimization/tms_opt_refined.py
   :language: python


MATLAB
''''''

.. literalinclude:: ../../simnibs/examples/optimization/tms_opt_refined.m
   :language: matlab


Auxiliary Dipole Method (ADM)
---------------------------------------

To use the Auxiliary Dipole Method (ADM), simply use a :file:`.ccd` or a :file:`.tcd` coil file that only contains dipole elements and set the :code:`method = 'ADM'`:

.. code-block:: python

  tms_opt.fnamecoil = os.path.join('legacy_and_other','Magstim_70mm_Fig8.ccd')
  tms_opt.method = 'ADM'

\

Python
''''''

.. literalinclude:: ../../simnibs/examples/optimization/tms_optimization_adm.py
   :language: python


MATLAB
''''''

.. literalinclude:: ../../simnibs/examples/optimization/tms_optimization_adm.m
   :language: matlab


Acknowledgements
------------------

We would like to thank Ole Numssen and Konstantin Weise for helping with the development of this
functionality, and Luis Gomez for contributing the code for the ADM optimization.

Further Reading
------------------
Please see :ref:`tmsoptimize_doc` for a detailed description of all TMS optimization options.

References
------------

`Weise, K., Numssen, O., Thielscher, A., Hartwigsen, G., Knösche, T.R. (in review) A novel approach to localize cortical TMS effects. bioRxiv, 595603. <https://doi.org/10.1101/595603>`_

`Gomez, L. J., Dannhauer, M., & Peterchev, A. V. (2020). Fast computational optimization of TMS coil placement for individualized electric field targeting. bioRxiv <https://doi.org/10.1101/2020.05.27.120022>`_



