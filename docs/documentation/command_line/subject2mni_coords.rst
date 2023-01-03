
.. _subject2mni_coords_docs:

subject2mni_coords
==================

Description
------------

Transforms coordinates from subject space to MNI space

Usage example
---------------

1. Suppose you open a simulation and see a hotspot at the position (-43, -17, 74). You want to determine where this hotspot is in MNI coordinates
2. Run

.. code-block:: bash

  subject2mni_coords -m m2m_ernie/ -c -43 -17 74

\
  This will output the corresponding coordinates in MNI space.

Further notes
--------------

* Type :code:`subject2mni_coords -h` for more information


