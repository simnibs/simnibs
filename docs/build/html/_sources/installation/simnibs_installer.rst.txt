.. _simnibs_installer:

Install SimNIBS
===============



A SimNIBS installation takes about 3 GB disk space.

.. note:: The installation procedure has changed in version 3.2. The previous simnibs installer app is now deprecated and will be removed in the near future.

.. note:: In case you have problems during the installation, please check :ref:`install_throubleshooting`

Windows
-------
1. `Download the SimNIBS installer <https://simnibs.drcmr.dk/userregistration2>`_

2. Double click the :file:`simnibs_installer_windows.exe` file. If a security warning shows up, click on *More info -> Run anyway*

3. Click through the installer wizard

4. Installation might take 5-10 minutes, depending on your computer and internet connection. Please be patient.
 

Testing the installation
'''''''''''''''''''''''''

Look for a shortcut called :code:`SimNIBS GUI` in your Start Menu


Linux
-----
1. `Download the SimNIBS installer <https://simnibs.drcmr.dk/userregistration2>`_


2. Run in a terminal 


  .. code-block:: bash
  
    cd ~/Downloads
    tar -xzf simnibs_installer_linux.tar.gz
    simnibs_installer/install

  \

3. Installation might take 5-10 minutes, depending on your computer and internet connection. Please be patient.

.. note:: The installer also has a silent mode (no GUI), type :code:`simnibs_installer/install -h` for more information


Testing the Installation
'''''''''''''''''''''''''
Start a new terminal window and type :code:`simnibs_gui`



MacOS
------
1. `Download the SimNIBS installer <https://simnibs.drcmr.dk/userregistration2>`_

2. Double click the :file:`simnibs_installer_macos.pkg` file.

3. Click through the installer wizard.

4. Installation might take 5-10 minutes, depending on your computer and internet connection. Please be patient.

.. note:: SimNIBS only supports MacOS versions ≥ 10.13 (High Sierra)


Testing the Installation
'''''''''''''''''''''''''
Open Launchpad and search for :code:`SimNIBS GUI`


Updating SimNIBS
-----------------

You can install the latest bugfix version of SimNIBS by starting a Command Prompt (Windows) or a Terminal (Linux and MacOS) window and typing

.. code-block:: bash

    update_simnibs

\

New feature versions require a new installation.

Software Dependencies
-----------------------
SimNIBS does not require any external dependencies for running simulations and post-processing operations on existing head models.

However, the head modelling pipelines have external dependencies.
  * :ref:`headreco_docs` requires MATLAB
  * :ref:`mri2mesh_docs` requires FSL and FreeSurfer

Please see :ref:`optional_deps` for more information on how to configure these dependencies
