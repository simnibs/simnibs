.. _simnibs_cli_docs:

simnibs
=========

Description
-------------

As alternative to start simulations. Simulations can be directly started from the command line using the *.mat* file with the simulation settings.

Usage example
--------------

1. Open a terminal and go to the directory of the “Ernie” example data set.
2. Run

.. code-block:: bash

  simnibs ernie_simu.mat

\

  You can run multiple TMS simulations at once by setting the :code:`--cpus` argument

.. code-block:: bash

  simnibs tms_simu.mat --cpus 3

\
  This will run 3 TMS simulations at once. However, each simulation might end up running slower.

