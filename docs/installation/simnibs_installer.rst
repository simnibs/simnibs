Using the SimNIBS Installer
===========================

As of version 3.0, SimNIBS is now distributed through the SimNIBS installer.


Windows
-------
1. `Download the SimNIBS installer <http://simnibs.drcmr.dk/>`_

2. Double click the :file:`install_simnibs` file

3. Select the desired version, accept the license and click *next*

4. Installation might take 10-20 minutes, depending on your computer and internet connection. Please be patient.
 
5. Log out and in again.

.. warning:: If you already have Anaconda or Miniconda installed in your computer, the links to the Anaconda Prompt and Anaconda Powershell in the Start Menu may be overwritten

Testing the installation
'''''''''''''''''''''''''

Look for a shortcut called :code:`SimNIBS GUI` in your Start Menu


Linux
-----
1. `Download the SimNIBS installer <http://simnibs.drcmr.dk/>`_

2. Extract the :file:`install_simnibs.tar.gz` file

3. Double click the :file:`install_simnibs` file

4. Select the desired version, accept the license and click on *next*

5. Installation might take 10-20 minutes, depending on your computer and internet connection. Please be patient.

.. warning:: The SimNIBS Installer does not support CentOS 6. However, it is possible to :ref:`install it using conda <conda-install>`

.. note:: The installer also has a silent mode, type :code:`./install_simnibs -h` for more information


Testing the Installation
'''''''''''''''''''''''''
Start a new terminal window and type :code:`simnibs_gui`


Updating SimNIBS
-----------------

If you already have SimNIBS >= 3.0 installed, it is possible to upgrade your
installation.

Just start a command prompt or a terminal and type :code:`update_simnibs`
