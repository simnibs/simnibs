.. _installation:

Installation
=============

Requirements
-------------

* SimNIBS requires a 64bit system. Approximately 8GB of memory are required to run smoothly.
 
* SimNIBS takes approximately 3GB disk space

* **Running simulations on existing head models LINK TO EXAMPLE DATASET does not require any external dependencies.**


Using the SimNIBS Installer
----------------------------


Windows
''''''''
SimNIBS has been tested on Windows 7 and Windows 10.

1. `Download the SimNIBS installer <http://simnibs.drcmr.dk/>`_

2. Double click the :file:`install_simnibs` file

3. Select the desired version, accept the license and click *next*
 
4. By default, SimNIBS is installed to :file:`C:\\Users\\<name>\\AppData\\Local\\SimNIBS`

5. Log out

6. You should now find shortcuts to SimNIBS in your Start Menu


Linux
'''''''

SimNIBS has been tested on Ubuntu 16.04 and 18.04, as well as CentOS 6 and 7.

1. `Download the SimNIBS installer <http://simnibs.drcmr.dk/>`_

2. Extract the :file:`install_simnibs.tar.gz` file

3. Double click the :file:`install_simnibs` file

4. Select the desired version, accept the license and click on *next*

5. By default, SimNIBS is installed to :file:`/home/<name>/SimNIBS`

6. Start a new terminal window

7. You should be able to start the GUI by typing :code:`simnibs_gui`


The installer also has a silent mode, type :code:`./install_simnibs -h` for more information


Using Conda (Advanced)
-----------------------

It is also possible to install SimNIBS using the `Conda <https://docs.conda.io/en/latest/>`_ package manager.


1. Go to out releases page `releases page <https://github.com/simnibs/simnibs/releases>`_ and find the version you want to install

2. Download the corresponding :file:`environment_win.yml` file (Windows) :file:`environment_unix.yml` file (Linux and OSX)

3. Create and activate the environment

  .. code-block:: bash
  
     conda env create -f environment_[win|unix].yml
     conda activate simnibs_env
  
  \

4. Install SimNIBS using :code:`pip`

  .. code-block:: bash
  
     pip install -f https://github.com/simnibs/simnibs/releases/tag/v<VERSION> simnibs

  \

5. (Optional) to install Menu icons and add SimNIBS to the system path, run the :code:`postinstall_simnibs` script. Type :code:`postinstall_simnibs -h` for help


Optional Dependencies
----------------------


MATLAB
'''''''
**Required by:** :ref:`headreco_docs` head segmentation pipeline

For :ref:`headreco_docs`, MATLAB needs to be configured such that typing :code:`matlab` on a terminal window will start it.

This is already the case in most Windows installations. However, on Linux and MacOSX you might need to create a link to the matlab executable somewhere in your system :code:`$PATH`

* Linux

  .. code-block:: bash
  
     sudo ln -s /usr/local/MATLAB/R<VERSION>/bin/matlab /usr/local/bin/matlab
  
  \

* MacOSX

  .. code-block:: bash
  
     sudo ln -s /Applications/MATLAB_R<VERSION>.app/bin/matlab /usr/local/bin/matlab
  
  \

If MATLAB is not installed in the default location, you can find out where it is installed by typing in a MATLAB terminal

.. code-block:: matlab

   matlabroot


SimNIBS also has a MATLAB API, available in the :file:`matlab/` subfolder of the SimNIBS installation directory.

FSL
'''
**Required by:** :ref:`mri2mesh_docs` head segmentation pipeline and :ref:`dwi2cond_docs` conductivity tensor reconstruction pipeline

1. Follow the instructions on `this link <http://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FslInstallation>`_

2. Add

.. code-block:: bash

   source /etc/fsl/5.0/fsl.sh

\
  in the end of the :file:`~/.bashrc` (Linux) or :file:`~/.bash_profile` (MacOSX) file (assuming that fsl is installed as usually into :file:`/etc/fsl/5.0`).


FreeSurfer
''''''''''
**Required by:** :ref:`mri2mesh_docs` head segmentation pipeline

1. Follow the instructions `here <http://freesurfer.net/fswiki/DownloadAndInstall>`_

2. Make sure that you have a registration file, and set the path in the :file:`~/.bashrc` (Linux) or :file:`~/.bash_profile` (MacOSX).



Uninstall SimNIBS
--------------------

Windows
'''''''
You can find the uninstaller in :code:`Add or Remove Programs`


Linux
'''''''
Run the :code:`uninstall simnibs`, located in the SimNIBS installation directory



Troubleshooting
-----------------

Please send an email to support@simnibs.org including the :file:`simnibs_install_log.txt` file, which can be found in your SimNIBS installation directory


