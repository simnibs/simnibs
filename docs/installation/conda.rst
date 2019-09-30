.. _conda-install:

Install Using Conda (Advanced)
===============================

It is also possible to install SimNIBS using the `Conda <https://docs.conda.io/en/latest/>`_ package manager.

1. Download and install the `Miniconda <https://docs.conda.io/en/latest/miniconda.html>`_ **Python 3** distribution.

2. Go to the latest `SimNIBS release page <https://github.com/simnibs/simnibs/releases/latest>`_ and download the :file:`environment_*.yml` file for your operating system

3. 
   * **Windows**: Open the the *Anaconda Prompt*, which can be found in the *Star Menu*.
   * **Linux** and **MacOS**: Open a terminal. Type

     .. code-block:: bash
  
        export PATH="$HOME/miniconda/bin:$PATH"
  
     \

     to ensure :code:`conda` is in your system PATH. You might need to change the command above depending on your miniconda install location.

4. Create and activate the :code:`simnibs_env` conda environment.

  .. code-block:: bash
  
     conda env create -f PATH/TO/environment_<OS>.yml
     conda activate simnibs_env
  
  \

  Where :file:`PATH/TO/environment_<OS>.yml` designates the path to the downloaded environment file, typically:

  * **Windows**:  :file:`%HOMEPATH%\\Downloads\\environment_win.yml` 
  * **Linux**:  :file:`~/Downloads/environment_linux.yml` 
  * **MacOS**:  :file:`~/Downloads/environment_macOS.yml` 

5. Install SimNIBS using :code:`pip`.

  .. code-block:: bash
  
     pip install -f https://github.com/simnibs/simnibs/releases/latest simnibs

  \

6. To setup the menu icons, file associations, the MATLAB library and add SimNIBS to the system path, run the :code:`postinstall_simnibs` script:

   * **Windows**

     .. code-block::
   
        md %LOCALAPPDATA%\SimNIBS
        postinstall_simnibs --copy-matlab --setup-links -d %LOCALAPPDATA%\SimNIBS

     \  

   * **Linux**

     .. code-block:: bash
   
        mkdir $HOME/SimNIBS
        postinstall_simnibs --copy-matlab --setup-links -d $HOME/SimNIBS

     \  

   * **MacOS**

     .. code-block:: bash
   
        mkdir -p $HOME/Applications/SimNIBS
        postinstall_simnibs --copy-matlab --setup-links -d $HOME/Applications/SimNIBS

     \  

   Type :code:`postinstall_simnibs -h` for more options.
