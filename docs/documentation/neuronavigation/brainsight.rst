.. _brainsight_doc:

Brainsight
==========
This module provides import and export functions for the `Brainsight <https://www.rogue-research.com/tms/brainsight-tms/>`_ TMS Navigator software.

Several coordinate systems are supported for export and import of coil data. SimNIBS always uses the `RAS` coordinate system (NIfTI). This and the `LPS` coordinate system used for DICOMS is supported here.

The Brainsight ecosystem provides two main ways of storing positions/orientations: `Targets` and `Samples`. Both are exported in a single `.txt` file.


How to use 
-----------

Export from Brainsight
######################
Choose `NIfTI:aligned` (preferred) or `World` coordinate systems while exporting the coil data.

Import to SimNIBS
#################
:code:`simnibs.brainsight.read(fn)` reads exported coil information from Brainsight and returns two :code:`simnibs.TMSLIST()` objects: One for `Targets` and one for `Samples`. The conversions from `LPS` (called `World` in Brainsight) to `RAS` for SimNIBS and the different coil axes definitions  are performed automatically.


..  code-block:: python
    :caption: Import a single Brainsight .txt file as :code:`simnibs.TMSLIST()`

    from simnibs import sim_struct, brainsight

    s = sim_struct.SESSION()

    fn = "exported_data.txt"
    tms_list_targets, tms_list_samples = brainsight().read(fn)  # read all Targets and Samples from file and return as TMSLIST() each
    s.add_tmslist(tms_list_targets)

    tms_list.pos[0].name  # <- name is filled with data from .txt.


Export from SimNIBS
###################
:code:`simnibs.brainsight.write(obj, fn)` writes an .xml file that is compatible for Brainsight import. The conversion between the different coil axes definitions is performed automatically.
When `World`, i.e. `LPS` coordinate system is selected for export the conversion is performed automatically.


.. code-block:: python
    :caption: Export a file for precomputed positions/orientations

    from simnibs import sim_struct, opt_struct, brainsight
    fn = "precomuted_coilpos.xml"

    ### export from TMSLIST
    tmslist = sim_struct.TMSLIST()
    tmslist.add_position()
    # ... define (multiple) positions ...
    brainsight().write(tmlist, fn, out_coord_space='LPS')

    ### export from POSITION
    pos = sim_struct.POSITOIN()
    pos.matsimnibs = ...
    brainsight().write(pos, fn) # out_coord_space default is 'RAS'

    ### export from np.ndarray / matsimnibs
    opt = opt_struct.TMSoptimize()
    # ... prepare optmization ...
    opt_mat = opt.run() # get optimal position
    brainsight().write(opt_mat, fn)

Notes
------

* The **same anatomical scan** has to be used for TMS Navigator and SimNIBS.
* The **same coil model** has to be used for field simulations and for real stimulation.
* Coordinate systems used to define coil axes for SimNIBS and Brainsight:

.. figure:: ../../images/coil_axesorientation_brainsight.png