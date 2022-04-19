.. _tdcsavoid_doc:

TDCSavoid
==========


Initialization
---------------

* **Python**

  .. code-block:: python

     from simnibs import optimization
     opt = optimization.TDCSoptimize()
     target = opt.add_avoid()



  \

* **MATLAB**

  .. code-block:: matlab

     opt = opt_struct('TDCSoptimize');
     opt.avoid(1)

  \ 


Attributes
-----------

* **positions**: *Nx3 list/array of floats (Python/MATLAB)*

  * **Desctiption**: Positions where the field is to be avoided (:ref:`see here for more
    information <positions_attribute_doc>`).


* **weight**: *float*

  * **Description**: Weight to be given to this avoid region. The largest the weight the
    more SimNIBS will avoid this region. Must be :math:`> 1`. 

  * **Default**: :math:`10^3`.

* **indexes**: *Nx1 list/array of ints (Python/MATLAB), optional*

  * **Description**: As an alternative to **positions**, you can select the **node**
    index or the **element** index (:ref:`see here for more
    information <indexes_attribute_doc>`).

  * **Default**: Get the points closest to the **positions**.

* **radius**: *float, optional*

  * :ref:`See here <radius_attribute_doc>`.

* **tissues**: *list/array of ints (Python/MATLAB), optional*

  * **Descrption**: List of tissue indices where the target is defined. Leave empty if
    all tissues in the leadfield.

  * **Default**: All tissues

  .. note:: You can leave **positions** and **indexes** empty and only assign **tissues**. In this case, the avoidance region will be the whole tissue. This can be useful, e.g. to avoid co-stimulation of the eyes when designing a tACS montage.



