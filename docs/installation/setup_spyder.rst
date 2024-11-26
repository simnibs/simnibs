.. _setup_spyder:

Setting SimNIBS up in Spyder (Advanced)
========================================

It is also possible to set up SimNIBS in the Spyder environment.

1. If you don't have Spyder installed yet, `download it <https://www.spyder-ide.org/download/>`_ and follow the installation instructions.

2. In Spyder *Tools-->Preferences* go to the *Python interpreter* tab to select the python interpreter of SimNIBS. 

.. image:: ../images/spyder_1.png
   :align: center
   :scale: 42 %
   
It is located in the *simnibs_env* subdirectory of the SimNIBS installation folder:

.. image:: ../images/spyder_2.png
   :align: center
   :scale: 60 %

3. Afterwards, start a new console in Spyder. This will now invoke the python interpreter of SimNIBS, and create an error message that states which packages need to be added to the SimNIBS installation:

.. image:: ../images/spyder_3.png
   :align: center
   :scale: 80 %
   
4. As last step, add these packages to the SimNIBS installation by running the following in a  :ref:`terminal window <windows_terminal>` (the stated version numbers may need to be updated):

  .. code-block::

	simnibs_python -m pip install spyder-kernels==3.0.*

  \

Starting again a new console in Spyder and importing SimNIBS should now work without errors:

  .. code-block::

	import simnibs
	print(simnibs.__version__)

  \
  