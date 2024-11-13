.. _tmslist_doc:

TMSLIST
=======

Describes a set of TMS simulations.


Initialization
---------------

* **Python**

  .. code-block:: python

     from simnibs import sim_struct
     s = sim_struct.SESSION()
     tms_list = s.add_tmslist()

  \

* **MATLAB**

  .. code-block:: matlab

     s = sim_struct('SESSION');
     s.poslist{1} = sim_struct('TMSLIST');

  \

Attributes
----------

.. _tmslist_fnamecoil:

* **fnamecoil**: *string (Python)/ character array (MATLAB)*

  * **Description**: Name of coil file. Coil files come in three types

    * *.nii.gz* files: NIfTI files with sampled magnetic vector potentials. Recommended, allows for faster simulations. (`Madsen et al., 2015 <https://doi.org/10.1016/j.brs.2015.07.035>`_)
    * *.ccd* files: Text files that describe the coil as a set of magnetic dipoles. Simulations with this type of coil are slower. (`Thielscher and Kammer, 2004 <https://doi.org/10.1016/j.clinph.2004.02.019>`_)
    * *.tcd* files: JSON based text files that describe the coil as a set of magnetic dipoles, line segments and/or sampled magnetic vector potentials.

  * **Examples**: *Python/MATLAB*

    Select the SimNIBS model for the Magstim 70mm figure-of-eight coil

    .. code-block:: matlab

       tmslist.fnamecoil = fullfile('legacy_and_other','Magstim_70mm_Fig8.ccd');

    \


  * **Note**: When using a :ref:`coil provided by SimNIBS <coil_fies>` you only need to use the file name. If using some other coil file, please use the full path.
  * **References**: `Madsen et al., 2015 <https://doi.org/10.1016/j.brs.2015.07.035>`_, `Thielscher and Kammer, 2004 <https://doi.org/10.1016/j.clinph.2004.02.019>`_

.. _tmslist_pos:

* **pos**: *list/array of POSITION structures (Python/MATLAB)*

  * **Description**: List of coil positions for the simulations
  * **Examples**: See the :ref:`documentation for the POSITION strucure <position_doc>`.

* **cond**: *list/array of COND structures (Python/MATLAB), optional*

  :ref:`Follow this link <cond_attribute_doc>`.

* **anisotropy_type**: *'scalar', 'vn', 'dir' or 'mc', optional*

  :ref:`Follow this link <anisotropy_type_attribute_doc>`.

* **aniso_maxratio**: *float*

  :ref:`Follow this link <aniso_maxratio_doc>`.

* **aniso_maxcond**: *float*

  :ref:`Follow this link <aniso_maxcond_doc>`.

* **solver_options**: *string (pytohn) / character array (MATLAB)*

  :ref:`Follow this link <solver_options_doc>`.




Examples
--------

* Set up a simulation with a coil over C3, pointing posteriorly.
  See the documentation on :ref:`session_doc` and the :ref:`position_doc` structures for more information.

  * *Python*

  .. code-block:: python

    from simnibs import sim_struct, run_simnibs
    # Create a SESSION structure
    S = sim_struct.SESSION()
    # Select the head mesh
    S.fnamehead = 'ernie.msh'
    # add a TMSLIST to the SESSION
    tms = S.add_tmslist()
    # Select the coil from those available in the coil_models subfolder
    tms.fnamecoil = os.path.join('legacy_and_other','Magstim_70mm_Fig8.ccd')
    # Add a new position
    pos = tms.add_position()
    # Place the coil over C3
    pos.centre = 'C3'
    # Point the coil towards CP3
    pos.pos_ydir = 'CP3'
    #  4 mm distance between coil and head
    pos.distance = 4

  \

  * *MATLAB*

  .. code-block:: matlab

    % Create a SESSION structure
    S = sim_struct('SESSION');
    % Select the head mesh
    S.fnamehead = 'ernie.msh';
    % Add a TMSLIST to the SESSION
    S.poslist{1} = sim_struct('TMSLIST');
    % Select the coil from those available in the -coil_models subfolder
    S.poslist{1}.fnamecoil = fullfile('legacy_and_other','Magstim_70mm_Fig8.ccd');
    % Place the coil over C3
    S.poslist{1}.pos(1).centre = 'C3';
    % Point the coil towards CP3
    S.poslist{1}.pos(1).pos_ydir = 'CP3';
    %  4 mm distance between coil and head
    S.poslist{1}.pos(1).distance = 4;

  \
