.. _convert_3_to_4:

convert_3_to_4
==============

Description
------------

Converts a head model created by *mri2mesh* or *headreco* for use in SimNIBS4

Usage example
-------------

1. Open a terminal and go to a directory with an old head model.
2. Run

.. code-block:: bash

  convert_3_to_4 m2m_<subID_old> m2m_<subID_new>

\

3. A folder m2m_<subID_new> will be created that contains the head model in SimNIBS4 convention

Further notes
---------------

* Type :code:`convert_3_to_4 -h` gives the help information
* When only the directory of the old head model is provided, the new folder will be named m2m_<subID_old>_v4
* Using the conversion function in python scripts is feasible by

.. code-block:: python

  from simnibs.cli.convert_3_to_4 import convert_old_new
  convert_old_new('m2m_<subID_old>', 'm2m_<subID_new>')

\