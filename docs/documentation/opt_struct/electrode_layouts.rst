.. _electrode_layouts:

Electrode and array layouts for TesFlexOptimization
================================================

TesFlexOptimization supports two layout types: *ElectrodeArrayPair* and *CircularArray*

ElectrodeArrayPair
-----------------------

Initialization
^^^^^^^^^^^^^^^

* **Python**

  .. code-block:: python

     from simnibs import opt_struct
     opt = opt_struct.TesFlexOptimization()
     electrode_layout = opt.add_electrode_layout("ElectrodeArrayPair")
     
  \

* **MATLAB**

  .. code-block:: matlab

     opt = opt_struct('TesFlexOptimization');
     opt.electrode{1} = opt_struct('ElectrodeArrayPair'); 

  \ 
  

Parameters and attributes
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autoclass:: simnibs.optimization.tes_flex_optimization.electrode_layout.ElectrodeArrayPair



CircularArray
-----------------------

Initialization
^^^^^^^^^^^^^^^

* **Python**

  .. code-block:: python

     from simnibs import opt_struct
     opt = opt_struct.TesFlexOptimization()
     electrode_layout = opt.add_electrode_layout("CircularArray")
     
  \

* **MATLAB**

  .. code-block:: matlab

     opt = opt_struct('TesFlexOptimization');
     opt.electrode{1} = opt_struct('CircularArray'); 

  \ 
  

Parameters and attributes
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autoclass:: simnibs.optimization.tes_flex_optimization.electrode_layout.CircularArray






