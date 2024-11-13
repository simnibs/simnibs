.. _regionofinterest_doc:

RegionOfInterest 
=================

Used for the definition of regions of interests on surfaces and in volumes. Please see also :file:`roi_definition.py` and :file:`roi_definition.m` in examples/utilities for usage examples.

Initialization
---------------

* **Python**

  .. code-block:: python

     from simnibs import opt_struct

     # TmsFlexOptimization
     opt = opt_struct.TmsFlexOptimization()
     roi = opt.add_region_of_interest()
     
     # TesFlexOptimization
     opt = opt_struct.TesFlexOptimization() 
     roi = opt.add_roi()

  \

* **MATLAB**

  .. code-block:: matlab

     % TmsFlexOptimization
     opt = opt_struct('TmsFlexOptimization');
     opt.roi{1} = opt_struct('RegionOfInterest');
     
     % TesFlexOptimization
     opt = opt_struct('TesFlexOptimization');
     opt.roi{1} = opt_struct('RegionOfInterest');
     
  \ 


Parameters and attributes
------------------------------------

.. autoclass:: simnibs.utils.region_of_interest.RegionOfInterest