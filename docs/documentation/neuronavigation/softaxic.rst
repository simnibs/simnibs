.. _softaxic_doc:

Softaxic
==========

How to use 
-----------

Import to SimNIBS
#################
:code:`simnibs.softaxic.read(fn)` reads exported coil information from SofTaxic and returns a :code:`simnibs.TMSLIST()` object.

..  code-block:: python
    :caption: Import a single SofTaxic .stmpx file as :code:`simnibs.TMSLIST()`

    from simnibs import sim_struct, softaxic

    s = sim_struct.SESSION()

    fn = "exported_data.stmpx"
    tms_list = softaxic().read(fn)  # read all targets from file and return as TMSLIST()
    s.add_tmslist(tms_list)

    tms_list.pos[0].name  # <- name is filled with pos-ID from .stmpx.


Notes
------

* export not available


