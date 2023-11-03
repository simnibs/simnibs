.. _get_fields_at_coordinates_doc:

get_fields_at_coordinates
=============================

Description
--------------

Interpolate field values at given points

Usage example
---------------

1. Suppose you want to accurately determine the electric field at a set of positions. The electric fields are in a mesh called *ernie_TDCS_1_scalar.msh*

2. Write the positions (x, y, z) where you want to sample the fields to a CSV file, as in :download:`this example <../../data/positions.csv>`

3. Run

.. code-block:: bash

  get_fields_at_coordinates -s positions.csv -m ernie_TDCS_1_scalar.msh

\

  This will output one file CSV for each field in the mesh, named for example *positions_E.csv*, *positions_magnE.csv*, …
  These files contain the fields (or field vectors) at each sampling position.

Further notes
-------------

* Type :code:`get_fields_at_coordinates -h` for more information and options
