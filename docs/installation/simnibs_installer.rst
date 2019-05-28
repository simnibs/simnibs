.. _simnibs_installer:

Install SimNIBS
===============

As of version 3.0, SimNIBS is now distributed through the SimNIBS installer.


Windows
-------
1. :ref:`Download the SimNIBS installer <download>`

2. Double click the :file:`install_simnibs_windows.exe` file

3. Select the desired version, accept the license and click *next*

4. Installation might take 15-30 minutes, depending on your computer and internet connection. Please be patient.
 
5. Log out and in again.

.. warning:: If you already have Anaconda or Miniconda installed in your computer, the links to the Anaconda Prompt and Anaconda Powershell in the Start Menu may be overwritten

Testing the installation
'''''''''''''''''''''''''

Look for a shortcut called :code:`SimNIBS GUI` in your Start Menu


Linux
-----
1. `Download the SimNIBS installer <http://simnibs.drcmr.dk/>`_

2. Run in a terminal 


  .. code-block:: bash
  
    cd ~/Downloads
    tar -xzf install_simnibs_linux.tar.gz
    install_simnibs/install_simnibs

  \

3. Select the desired version, accept the license and click on *next*

4. Installation might take 10-20 minutes, depending on your computer and internet connection. Please be patient.

.. warning:: The SimNIBS Installer does not support CentOS 6. However, it is possible to :ref:`install it using conda <conda-install>`

.. note:: The installer also has a silent mode, type :code:`./install_simnibs -h` for more information


Testing the Installation
'''''''''''''''''''''''''
Start a new terminal window and type :code:`simnibs_gui`



MacOS
------
1. `Download the SimNIBS installer <http://simnibs.drcmr.dk/>`_

2. Run in a terminal (use :code:`Cmd` + :code:`Space` and search for :code:`terminal`)


  .. code-block:: bash
  
    ~/Downloads/install_simnibs

  \

3. Select the desired version, accept the license and click on *next*

4. Installation might take 10-20 minutes, depending on your computer and internet connection. Please be patient.


Testing the Installation
'''''''''''''''''''''''''
Open Spotlight and search for :code:`SimNIBS`


Updating SimNIBS
-----------------

If you already have SimNIBS >= 3.0 installed, it is possible to upgrade your
installation.

Just start a Command Prompt (Windows) or a Terminal (Linux and MacOS) window and type

.. code-block:: bash

    update_simnibs

\
