.. _brainsight_doc:

Brainsight
==========
This module provides import and export functions for the `Brainsight <https://www.rogue-research.com/tms/brainsight-tms/>`_ TMS Navigator software.

SimNIBS requires an individual T1-NIfTI to be used during neuronavigation. This T1 scan should be the preprocessed T1 image that is created during the SimNIBS headmeshing procedure.
Although DICOMS can be used for the Brainsight neuronavigation, no import/export to SimNIBS is supported, due to (possible) different coordinate systems.

The Brainsight ecosystem provides two main ways of storing positions/orientations: `Targets` and `Samples`. Both are exported in a single `.txt` file. In addition, several coordinate systems can be chosen during export. Only `NIfTI:Aligned` is supported in the SimNIBS import.


How to use 
-----------

Exporting from Brainsight
######################
Choose `NIfTI:Aligned` coordinate systems while exporting the coil data.

Importing to SimNIBS
#################
:code:`simnibs.brainsight.read(fn)` reads exported coil position information files from Brainsight and returns two :code:`simnibs.TMSLIST()` objects: One for `Targets` and one for `Samples`. The conversion of the different coil axes definitions is performed automatically.


..  code-block:: python
    :caption: Import a single Brainsight .txt file as :code:`simnibs.TMSLIST()`

    from simnibs import sim_struct, brainsight

    s = sim_struct.SESSION()

    fn = "exported_data.txt"
    tms_list_targets, tms_list_samples = brainsight().read(fn)  # read all Targets and Samples from file and return as TMSLIST() each
    s.add_tmslist(tms_list_targets)

    tms_list.pos[0].name  # <- name is filled with data from .txt.

Exporting from SimNIBS
###################
:code:`simnibs.brainsight.write(obj, fn)` writes a .txt file that is compatible for Brainsight import. The conversion between the different coil axes definitions is performed automatically.


.. code-block:: python
    :caption: Export a file for precomputed positions/orientations

    from simnibs import sim_struct, opt_struct, brainsight
    fn = "precomuted_coilpos.xml"

    ### export from TMSLIST
    tmslist = sim_struct.TMSLIST()
    tmslist.add_position()
    # ... define (multiple) positions ...
    brainsight().write(tmlist, fn)

    ### export from POSITION
    pos = sim_struct.POSITION()
    pos.matsimnibs = ...
    brainsight().write(pos, fn)

    ### export from np.ndarray / matsimnibs
    opt = opt_struct.TMSoptimize()
    # ... prepare optmization ...
    opt_mat = opt.run() # get optimal position
    brainsight().write(opt_mat, fn)

Notes
-----

* The **same anatomical scan** has to be used for Brainsight and SimNIBS.
* The **same coil model** has to be used for field simulations and for real stimulation.
* Coordinate systems used to define coil axes for SimNIBS and Brainsight:

.. figure:: ../../images/coil_axesorientation_brainsight.png