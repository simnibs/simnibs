.. _simnibs_installer:

Install SimNIBS
===============



A SimNIBS installation takes about 3 GB disk space.

.. note:: In case you have problems during the installation, please check :ref:`install_throubleshooting`

Windows
-------
1. `Download the SimNIBS installer <https://github.com/simnibs/simnibs/releases/latest>`_

2. Double click the :file:`simnibs_installer_windows.exe` file. If a security warning shows up, click on *More info -> Run anyway*

3. Click through the installer wizard

4. Installation might take 5-10 minutes, depending on your computer and internet connection. Please be patient.
 

Testing the installation
'''''''''''''''''''''''''

Look for a shortcut called :code:`SimNIBS GUI` in your Start Menu


Linux
-----
1. `Download the SimNIBS installer <https://github.com/simnibs/simnibs/releases/latest>`_


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

.. note:: The SimNIBS GUI fails on some linux distributions, in particular those with Wayland. Workaround: ld preloading of libstdc++.so.6 seems to help. Example: export LD_PRELOAD=/usr/lib/libstdc++.so.6 (path needs to be adjusted according to library path on local system)


MacOS
------
1. `Download the SimNIBS installer <https://github.com/simnibs/simnibs/releases/latest>`_

2. Double click the :file:`simnibs_installer_macos.pkg` file.

3. Click through the installer wizard.

4. Installation might take 5-10 minutes, depending on your computer and internet connection. Please be patient.

.. note:: SimNIBS only supports Apple Silicon (Support for Intel-based Macs has been deprecated starting from version 4.5)


Testing the Installation
'''''''''''''''''''''''''
Open Launchpad and search for :code:`SimNIBS GUI`


Software Dependencies
-----------------------
SimNIBS does not require any external dependencies for creating head models, running simulations and post-processing operations. The preparation of conductivity tensors for GM and WM from diffusion MRI data using :ref:`dwi2cond_docs` requires FSL. Please see :ref:`optional_deps` for more information.
