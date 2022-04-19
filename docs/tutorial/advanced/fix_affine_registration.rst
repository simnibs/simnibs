.. _fix_affine_registration_tutorial:

.. role:: python(code)
   :language: python

Fixing the Affine Registration
==============================

.. note:: This tutorial applies only when setting :python:`init_type = "atlas"` in the :ref:`charm_docs` configuration file (settings.ini).

When encountering spurious segmentation results, we have found that this can often be traced back to a suboptimal affine registration between the structural scan(s) and the atlas. Identifying what exactly is causing the segmentation to fail can be a bit tricky (at least initially) so here we provide a few examples showing the type of errors that one might encounter and how they can be alleviated by modifying the settings of :ref:`charm_docs`. Charm writes a `settings.ini` inside the m2m-folder of the subject. Make a copy of this and then modify it as described below. This way, different changes can be made for individual subjects. To use the modified settings file, pass it using the `--usesettings` argument when running :ref:`charm_docs`.

Here we use the visualizations provided by :ref:`charm_docs` in the `charm_report.html` document, specifically, the `Coregistered template` and `Tissue labels` overlays.

1. Bad Rotation
---------------

In the final segmentation, we see what appears to be a downwards rotation around the *x*-axis. We can fix this by changing the default value of :python:`affine_rotations` (which is :python:`[-7, -3.5, 0, 3.5, 7]`) to for example :python:`affine_rotations = [0]`, thus restricting the rotational initializations around the *x*-axis.

:python:`affine_rotations` are specified as a list of initial rotations in degrees.

.. figure:: ../../images/tutorial_fix_affine_registration/rotation_issue_affine.png

       Fig 1.1. Affine (rotation issue).

.. figure:: ../../images/tutorial_fix_affine_registration/rotation_issue_final.png

       Fig 1.2. Final segmentation (rotation issue).

.. figure:: ../../images/tutorial_fix_affine_registration/rotation_issue_affine_fixed.png

       Fig 1.3. Affine (rotation fixed).

.. figure:: ../../images/tutorial_fix_affine_registration/rotation_issue_final_fixed.png

       Fig 1.4. Final segmentation (rotation fixed).


1. Bad Scaling
--------------

In the final segmentation, we see a clear dent in the frontal part of the head. Looking at the affine registration it seems that the atlas was scaled down too much. We can fix this by changing the default value of :python:`affine_scales` (which is :python:`[[0.85, 0.85, 0.85], [0.9, 0.9, 0.85], [0.95, 0.95, 0.85]]`) to for example :python:`affine_scales = [[0.9, 1.0, 0.85], [0.95, 1.0, 0.85]]` where we make the initializations of the scaling (along the *x*- and *y*-axis) larger.

:python:`affine_scales` is specified as a list of initial scalings in along each axis, i.e., [x, y, z].

.. figure:: ../../images/tutorial_fix_affine_registration/scale_issue_affine.png

       Fig 2.1. Affine (scale issue).

.. figure:: ../../images/tutorial_fix_affine_registration/scale_issue_final.png

       Fig 2.2. Final segmentation (scaling issue).

.. figure:: ../../images/tutorial_fix_affine_registration/scale_issue_affine_fixed.png

       Fig 2.3. Affine (scaling fixed).

.. figure:: ../../images/tutorial_fix_affine_registration/scale_issue_final_fixed.png

       Fig 2.4. Final segmentation (scaling fixed).


3. Bad Neck Deformation
-----------------------

In the final segmentation, we see that the neck is deformed too far back. This can be fixed by restricting the posterior search bound of the neck deformation. Specifically, we can replace the default value of :python:`neck_search_bounds` (which is :python:`[-0.3, 0.1]`) with :python:`neck_search_bounds = [0, 0.2]` to prevent posterior deformation and allow slightly more deformation in the anterior direction.

:python:`neck_search_bounds` is specified as a list with a posterior and an anterior search bound.

.. figure:: ../../images/tutorial_fix_affine_registration/neck_issue_affine.png

       Fig 3.1. Affine (neck deformation issue).

.. figure:: ../../images/tutorial_fix_affine_registration/neck_issue_final.png

       Fig 3.2. Final segmentation (neck deformation issue).

.. figure:: ../../images/tutorial_fix_affine_registration/neck_issue_affine_fixed.png

       Fig 3.3. Affine (neck deformation fixed). Since the neck deformation is performed *after* the affine registration, there is no visual effect of this,  but the result can be appreciated in the final segmentation (fig. 3.4).

.. figure:: ../../images/tutorial_fix_affine_registration/neck_issue_final_fixed.png

       Fig 3.4. Final segmentation (neck deformation fixed)
