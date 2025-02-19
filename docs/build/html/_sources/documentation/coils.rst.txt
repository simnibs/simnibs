.. _coil_fies:


Coil Models Included in SimNIBS
==================================


Default Coil Files
-------------------

Please use the files in the coil_models subfolder *Drakaki_BrainStim_2022*. They were reconstructed from direct measurements of their fields, making them quite accurate. Further details, including maximal dI/dt values at 100% MSO, are listed in the corresponding publication.

Reference
''''''''''

`Drakaki M, Mathiesen C, Siebner HR, Madsen KH, Thielscher A. (2022). Database of 25 validated coil models for electric field simulations for TMS. Brain Stimulation. <https://doi.org/10.1016/j.brs.2022.04.017>`_

Flexible Coil Files
-------------------

Models of the flexible H1, H4 and H7 coils from Brainsway and the MST-Twin coil with movable parts from MagVenture can be found in subfolder *flexible_coils*. The Brainsway coils were constructed based on X-ray images.
A description of the TMS coil definition file format (.tcd) can be found in :ref:`file_formats`.

Reference
''''''''''

`Worbs T, Rumi B, Madsen KH, Thielscher A, Personalized electric field simulations of deformable large TMS coils based on automatic position and shape optimization, bioRxiv <https://doi.org/abc>`_

Previous Coil Files
-------------------

The two previously included models of the Magstim 70mm and MagVenture MC-B70 Figure-of-Eight coils can be found in subfolder *legacy_and_other*. They were constructed based on X-ray images.


Reference
''''''''''

`Thielscher, Axel, and Thomas Kammer. "Electric field properties of two commercial figure-8 coils in TMS: calculation of focality and efficiency." Clinical neurophysiology 115.7 (2004): 1697-1708. <https://doi.org/10.1016/j.clinph.2004.02.019>`_


Extra Coil Files
----------------

Please see :download:`here <Deng_Brain_Stim_2013_docu.pdf>` for more information. These coil files can be downloaded after installation with the command

.. code-block::

  download_coils

We would like to thank Zhi-De Deng for sharing his coil definition files with us and Mursel Karadas for converting them to SimNIBS format.
These models were constructed from geometrical values about the coil windings taken from literature.
Please note that some of the models represent non-planar coils.
Setting the coil position for those in the SimNIBS GUI can easily result in the simulation of physically impossible coil placements,
i.e. corresponding to parts of the coil being inside the head, and are only meant for expert use!

Note (29.8.22)
''''''''''''''

The coil model "No9_H1.nii.gz" is inaccurate, mainly (but not only) because it is left-right mirrored. We are looking into ways to provide a more accurate model of that coil.


Reference
''''''''''

`Deng, Zhi-De, Sarah H. Lisanby, and Angel V. Peterchev. "Electric field depth–focality tradeoff in transcranial magnetic stimulation: simulation comparison of 50 coil designs." Brain stimulation 6.1 (2013): 1-13. <https://doi.org/10.1016/j.brs.2012.02.005>`_ 
