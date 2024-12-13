.. _scripting_tutorial:

Scripting Simulations
======================

This tutorial is based on the `Example Dataset <https://github.com/simnibs/example-dataset/releases/latest/download/simnibs4_examples.zip>`_.
Please download it before continuing.

SimNIBS offers both *Python* and *MATLAB* interfaces for setting up and running simulations.
In both, we use a set of nested structures to define the simulation, and then use the *run_simnibs* function to run the simulation.

.. _run_scripts:

Running Scripts
-----------------

Python
''''''

Typing 
  .. code-block:: bash

   simnibs_jupyter

\  
in a terminal window will start up JupyterLab with the correct *Python* interpreter. Alternatively, you can use :ref:`simnibs_python <simnibs_python_cli>` to run SimNIBS *Python* scripts. In case you have Spyder installed, you can also :ref:`set it up to use the SimNIBS installation <setup_spyder>`.

MATLAB
''''''

Add the SimNIBS *MATLAB* functions to the *MATLAB* path. In default installations, you can call

* **Windows:**

  .. code-block:: matlab

    addpath('C:\Users\<USER_NAME>\AppData\Local\SimNIBS\matlab_tools')

\

* **Linux:**

  .. code-block:: matlab

    addpath('/home/<USER_NAME>/SimNIBS/matlab_tools')

\

* **MacOS**

  .. code-block:: matlab

    addpath('/Users/<USER_NAME>/Applications/SimNIBS/matlab_tools')

\


Starting a SESSION and Selecting a Head Mesh
---------------------------------------------

The base structure for SimNIBS scripts is the :ref:`session_doc`. It may contain many
simulations of different types (TMS or tDCS), sharing the same head model.

We always start our scripts by initializing a :ref:`session_doc` *class* (*Python*) or *struct* (*MATLAB*), selecting a head mesh and the output folder.
Here, we will assume that the scripts are placed in the same directory as the :file:`m2m_ernie` directory of the example dataset.
If the scripts are not in the same folder as the subject folder, you should also give the path to the subject folder.


* *Python*

.. code-block:: python

    from simnibs import sim_struct, run_simnibs

    # Initalize a session
    s = sim_struct.SESSION()
    # Name of head mesh
    s.subpath = 'm2m_ernie'
    # Output folder
    s.pathfem = 'tutorial/'

* *MATLAB*

  .. code-block:: matlab

    % Initialize a session
    s = sim_struct('SESSION');
    % Name of head mesh
    s.subpath = 'm2m_ernie';
    % Output folder
    s.pathfem = 'tutorial/';


.. seealso:: Output and post-processing options are also configured in the :ref:`session_doc` structure. Please see the :ref:`documentation <session_doc>` for more details.


Setting up a TMS Simulation
----------------------------


Now, we want to set-up a TMS simulation.
To do it, we add a :ref:`tmslist_doc` to the :ref:`session_doc` structure and select a coil model (:ref:`list of available coils <coil_fies>`).


* *Python*

  .. code-block:: python

     # Initialize a list of TMS simulations
     tmslist = s.add_tmslist()
     # Select coil
     tmslist.fnamecoil = os.path.join('legacy_and_other','Magstim_70mm_Fig8.ccd')


* *MATLAB*

  .. code-block:: matlab

    % Initialize a list of TMS simulations
    s.poslist{1} = sim_struct('TMSLIST');
    % Select coil
    s.poslist{1}.fnamecoil = fullfile('legacy_and_other','Magstim_70mm_Fig8.ccd');

Now we need to set a position for our coil. Suppose we want to place it on position C1, pointing
posteriorly. You can do it by

* *Python*

  .. code-block:: python

     # Initialize a coil position
     pos = tmslist.add_position()
     # Select coil centre
     pos.centre = 'C1'
     # Select coil direction
     pos.pos_ydir = 'CP1'


* *MATLAB*

  .. code-block:: matlab

    % Select coil centre
    s.poslist{1}.pos(1).centre = 'C1';
    % Select coil direction
    s.poslist{1}.pos(1).pos_ydir = 'CP1';


We can set many coil positions to a single :ref:`tmslist_doc`. For example, we can add one
more coil position, now with the coil pointing towards Cz.


* *Python*

  .. code-block:: python

     # Add another position
     pos_superior = tmslist.add_position()
     # Centred at C1
     pos_superior.centre = 'C1'
     # Pointing towards Cz
     pos_superior.pos_ydir = 'Cz'


* *MATLAB*

  .. code-block:: matlab

    % Centred at C1
    s.poslist{1}.pos(2).centre = 'C1';
    % Pointing towards Cz
    s.poslist{1}.pos(2).pos_ydir = 'Cz';



.. seealso:: Coil positions are set through the  :ref:`position_doc` structure. It also allows you to set stimulator intensity (dI/dt values) and define coil positions in other ways. Please see the :ref:`documentation <position_doc>` for more information.


Setting up a tDCS Simulation
-----------------------------

To perform a tDCS simulation, we begin by setting a :ref:`tdcslist_doc` structure to the :ref:`session_doc` and setting the current flow through each channel. Here, we will only use two electrodes and set the current to 1mA. The first electrode will be a cathode, and the second an anode.

* *Python*

  .. code-block:: python

     # Initialize a tDCS simulation
     tdcslist = s.add_tdcslist()
     # Set currents
     tdcslist.currents = [-1e-3, 1e-3]


* *MATLAB*

  .. code-block:: matlab

    % Initialize a tDCS simulation
    s.poslist{2} = sim_struct('TDCSLIST');
    % Set currents
    s.poslist{2}.currents = [-1e-3 1e-3];

Let's first set the cathode. Suppose we want a 70x50mm rectangular over C3, pointing towards Cz.


* *Python*

  .. code-block:: python

     # Initialize the cathode
     cathode = tdcslist.add_electrode()
     # Connect electrode to first channel (-1e-3 mA, cathode)
     cathode.channelnr = 1
     # Electrode dimension
     cathode.dimensions = [50, 70]
     # Rectangular shape
     cathode.shape = 'rect'
     # 5mm thickness
     cathode.thickness = 5
     # Electrode Position
     cathode.centre = 'C3'
     # Electrode direction
     cathode.pos_ydir = 'Cz'


* *MATLAB*

  .. code-block:: matlab

     % Connect electrode to first channel (-1e-3 mA, cathode)
     s.poslist{2}.electrode(1).channelnr = 1;
     % Electrode dimension
     s.poslist{2}.electrode(1).dimensions = [50 70];
     % Rectangular shape
     s.poslist{2}.electrode(1).shape = 'rect';
     % 5mm thickness
     s.poslist{2}.electrode(1).thickness = 5;
     % Electrode Position
     s.poslist{2}.electrode(1).centre = 'C3';
     % Electrode direction
     s.poslist{2}.electrode(1).pos_ydir = 'Cz';


Now we need to configure the anode. Let's set a 30x30mm circular electrode over C4.

* *Python*

  .. code-block:: python

     # Add another electrode
     anode = tdcslist.add_electrode()
     # Assign it to the second channel
     anode.channelnr = 2
     # Electrode diameter
     anode.dimensions = [30, 30]
     # Electrode shape
     anode.shape = 'ellipse'
     # 5mm thickness
     anode.thickness = 5
     # Electrode position
     anode.centre = 'C4'


* *MATLAB*

  .. code-block:: matlab

     % Assign the electrode to the second channel
     s.poslist{2}.electrode(2).channelnr = 2;
     % Electrode diameter
     s.poslist{2}.electrode(2).dimensions = [30 30];
     % Electrode shape
     s.poslist{2}.electrode(2).shape = 'ellipse';
     % Electrode thickness
     s.poslist{2}.electrode(2).thickness = 5;
     % Electrode position
     s.poslist{2}.electrode(2).centre = 'C4';


.. seealso:: Electrodes are defined through the highly flexible :ref:`electrode_struct_doc` structure. Please see the :ref:`documentation <electrode_struct_doc>` for more information. Please note that it is also possible to connect multiple electrodes to a single channel, which is not possible to do in the GUI.

Running Simulations
---------------------

After the simulations are set, we can use the *run_simnibs* function to run the
simulations:

.. code-block:: matlab

   run_simnibs(s)


Now run the script in *Python* (using the :ref:`simnibs_python <simnibs_python_cli>` command) or in *MATLAB*.
After the simulations have finished running, the results can be found in the newly created
:file:`tutorial/` folder.

* Download the full :download:`Python <../data/tutorial_python.py>` and :download:`MATLAB <../data/tutorial_matlab.m>` scripts.


More Examples
----------------

More examples can be found in the :file:`examples/` folder in your SimNIBS installation directory. In default installations, it can be found at

* **Windows:**

  :file:`C:\\Users\\<USER_NAME>\\AppData\\Local\\SimNIBS\\examples`

* **Linux:**

  :file:`/home/<USER_NAME>/SimNIBS/examples`

* **MacOS:**

  :file:`/Users/<USER_NAME>/Applications/SimNIBS.app/examples`

Further Reading
----------------

* Tutorial on :ref:`visualization_tutorial`
* More information on the :ref:`sim_struct_doc`
* For an example on how to do group analysis in SimNIBS, please see the `SimNIBS 2.1 tutorial paper <https://doi.org/10.1101/500314>`_.
