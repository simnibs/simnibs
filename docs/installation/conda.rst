.. _conda-install:

Install Using Conda (Advanced)
===============================

It is also possible to install SimNIBS using the `Conda <https://docs.conda.io/en/latest/>`_ package manager.

.. note::
  If the pip installation fails, please locate the appropriate wheel for your system on the release page (e.g., the `latest <https://github.com/simnibs/simnibs/releases/latest>`_), copy its address, and feed this directly to ``pip install``, e.g.,

    pip install https://github.com/simnibs/simnibs/releases/download/v4.0.0/simnibs-4.0.0-cp39-cp39-win_amd64.whl

Windows
--------

1. Download and install the `Miniconda <https://docs.conda.io/en/latest/miniconda.html>`_ **Python 3** distribution.

2. Download the `SimNIBS Windows environment file <https://github.com/simnibs/simnibs/releases/latest/download/environment_win.yml>`_

3. Open  *Anaconda Prompt*, which can be found in the *Start Menu*.

4. Run in the Prompt:

  .. code-block:: bash

      conda env create -f "%USERPROFILE%\Download\environment_win.yml"
      conda activate simnibs_env
      pip install -f https://github.com/simnibs/simnibs/releases/latest simnibs

  \

5. (Optional) To setup the menu icons, file associations, the MATLAB library and add SimNIBS to the system path, run the :code:`postinstall_simnibs` script:

  .. code-block::

     md "%USERPROFILE%\SimNIBS"
     postinstall_simnibs --setup-links -d "%USERPROFILE%\SimNIBS"

  \

Linux
-------

1. Download and install the `Miniconda <https://docs.conda.io/en/latest/miniconda.html>`_ **Python 3** distribution.

2. Download the `SimNIBS Linux environment file <https://github.com/simnibs/simnibs/releases/latest/download/environment_linux.yml>`_

3. Run in a terminal window:

  .. code-block:: bash

      export PATH="$HOME/miniconda/bin:$PATH" # This part can change depending on your miniconda installation
      conda env create -f ~/Downloads/environment_linux.yml
      conda activate simnibs_env
      pip install -f https://github.com/simnibs/simnibs/releases/latest simnibs

  \

4. (Optional) To setup the menu icons, file associations, the MATLAB library and add SimNIBS to the system path, run the :code:`postinstall_simnibs` script:

  .. code-block:: bash

     mkdir $HOME/SimNIBS
     postinstall_simnibs --setup-links -d $HOME/SimNIBS

  \


MacOS
------

1. Download and install the `Miniconda <https://docs.conda.io/en/latest/miniconda.html>`_ **Python 3** distribution.

2. Download the `SimNIBS OSX environment file <https://github.com/simnibs/simnibs/releases/latest/download/environment_macOS.yml>`_

3. If you are using an Intel-based Mac (x86_64), run the following in a terminal window:

  .. code-block:: bash
  
      export PATH="$HOME/miniconda/bin:$PATH" # This part can change depending on your miniconda installation
      conda env create -f ~/Downloads/environment_macOS.yml
      conda activate simnibs_env
      pip install -f https://github.com/simnibs/simnibs/releases/latest simnibs
  
  \

If you are using an ARM-based Mac, you will have to create an x86_64 environment in order to build simnibs. Run this instead

  .. code-block:: bash
  
      export PATH="$HOME/miniconda/bin:$PATH" # This part can change depending on your miniconda installation
      conda env create -f ~/Downloads/environment_macOS.yml --platform osx-64
      conda activate simnibs_env
      pip install -f https://github.com/simnibs/simnibs/releases/latest simnibs
  
  \

If the :code:`--platform` argument is not available in your version of conda, use the following workaround
  
  .. code-block:: bash
  
      export PATH="$HOME/miniconda/bin:$PATH" # This part can change depending on your miniconda installation
      CONDA_SUBDIR=osx-64 conda env create -f ~/Downloads/environment_macOS.yml
      conda env config vars set CONDA_SUBDIR=osx-64 -n simnibs_env
      conda activate simnibs_env
      pip install -f https://github.com/simnibs/simnibs/releases/latest simnibs
  
  \

4. (Optional) To setup the menu icons, file associations, the MATLAB library and add SimNIBS to the system path, run the :code:`postinstall_simnibs` script:

  .. code-block:: bash

     mkdir -p $HOME/Applications/SimNIBS
     postinstall_simnibs --setup-links -d $HOME/Applications/SimNIBS

  \
