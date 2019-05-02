.. _installation:

Installation
=============

Overview
---------

* `Download SimNIBS <http://simnibs.drcmr.dk/>`_

* Download the example data set **LINK TO EXAMPLE**. 

* Installation requires an internet connection.

* SimNIBS now has two programs to create individual head models from MRI data

  * :ref:`mri2mesh_docs` is an older program. It is not available on Windows and requires
    FreeSurfer (version 5.3.0 or newer) and FSL (version 5.0.9 or newer).
  * :ref:`headreco_docs` is the new head modelling pipeline. It is available for all major
    operating systems and requires MATLAB to run.

* SimNIBS 3.0 has been tested in Ubuntu 18.04, Ubuntu 16.04, CentOS 6, CentOS 7, MacOS 10.13, Windows 7 and Windows 10. A 64-bit computer and 8 GB memory are recommended.

* In case of problems during the installation, contact support@simnibs.org


* SimNIBS is copyrighted |copy| by its :ref:`authors <contributors>` and licensed under :download:`GPL v3 <../LICENSE.txt>`.

.. |copy|   unicode:: U+000A9 .. COPYRIGHT SIGN

Windows
--------
1. `Download SimNIBS <http://simnibs.drcmr.dk/>`_
2. Double click in the auto-installer.
3. SimNIBS will unpack the files and afterwards start a PowerShell window. **Do not
   click on this window**. This step can take between 10 and 30 minutes, depending on
   your computer and internet connection.

Uninstall
'''''''''''

To uninstall SimNIBS on Windows, double click on the **Uninstall SimNIBS** in the SimNIBS directory

Linux and MacOSX
------------------

1. `Download SimNIBS <http://simnibs.drcmr.dk/>`_
2. Start a terminal window and navigate to the download folder
3. Uncompress the SimNIBS folder

.. code-block:: bash

      tar -zxvf simnibs-3.*.tar.gz

\
4. Go to simnibs-3.X.X-Linux64 folder

.. code-block:: bash

      cd simnibs-3.*

\
5. Execute the automated installation script

.. code-block:: bash

      ./install_simnibs.sh

\
6. Open a new terminal window or tab, a message should appear

.. code-block:: text

      Setting up paths for SimNIBS 3.X.X
      SIMNIBSDIR /path/to/simnibs

\

7. (MacOSX only) Download and install `Gmsh <http://gmsh.info/>`_.

.. note:: The installation script will also modify the ~/.bashrc file (on Linux) or the ~/.bash_profile file (on MacOSX) to set-up the SimNIBS environment whenever a new terminal window starts.


.. _matlab_headreco:

(Optional) Configuring MATLAB for **headreco**
'''''''''''''''''''''''''''''''''''''''''''''''''

In order to use the **headreco** head meshing script, MATLAB needs to be configured such that typing **matlab** on a terminal window will start it.
To configure it, type in a terminal window

On Linux:

.. code-block:: text

   sudo ln -s /usr/local/MATLAB/R<VERSION>/bin/matlab /usr/local/bin/matlab

On MacOSX:

.. code-block:: text

   sudo ln -s /Applications/MATLAB_R<VERSION>.app/bin/matlab /usr/local/bin/matlab

Substitute <VERSION> with your MATLAB version.

If MATLAB is not installed in the default location, you can find out where it is installed by typing in a MATLAB terminal

.. code-block:: text

   matlabroot

(Optional) installing FSL for **mri2mesh** and **dwi2cond**
'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

FSL is needed to run the head segmentation script :ref:`mri2mesh_docs` and the DWI to conductivity tensors script :ref:`dwi2cond_docs`

1. Follow the instructions on `this link <http://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FslInstallation>`_

2. Add

.. code-block:: text

   source /etc/fsl/5.0/fsl.sh

\
  in the end of the :file:`~/.bashrc` (Linux) or :file:`~/.bash_profile` (MacOSX) file (assuming that fsl is installed as usually into :file:`/etc/fsl/5.0`).


(Optional) installing FreeSurfer for **mri2mesh**
''''''''''''''''''''''''''''''''''''''''''''''''''

FreeSurfer is needed to run the head segmentation script :ref:`mri2mesh_docs`

1. Follow the instructions `here <http://freesurfer.net/fswiki/DownloadAndInstall>`_

2. Make sure that you have a registration file, and set the path in the :file:`~/.bashrc` (Linux) or :file:`~/.bash_profile` (MacOSX).




Troubleshooting
''''''''''''''''

* If the message

.. code-block:: text

      Setting up paths for SimNIBS 3.X.X
      SIMNIBSDIR /path/to/simnibs

\

  does not show up when starting a new terminal, please ensure you are using a **bash** terminal (and **not tcsh**)  

* If you are having problems installing with **sudo**, try adding the **-E** option

.. code-block:: bash

      sudo -E ./install_simnibs.sh

\

* When installing on older Linux distribution such as **CentOS 6**, you might need to
  download an `older Gmsh binary <gmsh.info/bin/Linux/gmsh-2.16.0-Linux64.tgz>`_ and
  substitute in the :file:`$SIMNIBSDIR/bin` folder.

Uninstall
'''''''''

Linux
^^^^^^
.. code-block:: bash

  rm -r $SIMNIBSDIR
  sed -i.bak '/SIMNIBS/d' ~/.bashrc

\


MacOSX
^^^^^^^
.. code-block:: bash

  rm -r $SIMNIBSDIR
  sed -i.bak '/SIMNIBS/d' ~/.bash_profile

\


(Optional) MATLAB Functions
-----------------------------

MATLAB functions for interacting with SimNIBS are available in the :file:`matlab/` subfolder of the SimNIBS installation directory.

