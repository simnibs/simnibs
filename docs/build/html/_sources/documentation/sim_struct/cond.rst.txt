.. _cond_struct_doc:

COND
====

Describes tissue conductivity values

Initialization
--------------

* **Python**

  .. code-block:: python

     from simnibs import sim_struct
     s = sim_struct.SESSION()
     tdcs_list = s.add_tdcslist()
     gm_cond = tdcs_list.cond[1]

  \

* **MATLAB**

  .. code-block:: matlab

     s = sim_struct('SESSION');
     s.poslist{1} = sim_struct('TDCSLIST');
     gm_cond = s.poslist{1}.cond(2)

  \ 

Attributes
----------

* **value**: *float*

  * **Description**: Conductivity value, in S/m

* **Name**: *string (Python)/character array (MATLAB), optional*


  * **Description**: Name of the tissue

* **distribution_type**: *'uniform', 'normal', 'beta' or None/unset (Python/MATLAB), optional*
   * **Description**: type of distribution for Uncertainty Quantification. Default: *None/unset*

   .. warning:: Setting this property will trigger a run of the UQ algorithm. Please see the tutorial on :ref:`uq_tutorial` for more information

\

* **distribution_parameters**: *list/array of floats (Python/MATLAB)*
   * **Description**: Sets the parameters for distribution defined in **distribution_type**.

     * if distribution_type is *'uniform'*: [min_value, max_value]
     * if distribution_type is *'normal'*: [mean, standard_deviation]
     * if distribution_type is *'beta'*: [p, q, min_value, max_value]


Examples
----------

.. note:: In this examples, we will consider a list/array of COND structures

* Change the value of gray matter conductivity. See the :ref:`Standard conductivity values <conductivity>`.

  * *Python*
    
    .. code-block:: python

       # In Python, we have to take the cond struct
       # with the index tissue_number - 1
       cond[1].value = 0.2

    \

  * *MATLAB*

    .. code-block:: matlab

       cond(2).value = 0.2;

    \

* Add the conductivity of 10 S/m and a name for a new tissue, labeled as 11

  * *Python*
    
    .. code-block:: python

       cond[10].name = 'Awesome tissue'
       cond[10].value = 10

    \

  * *MATLAB*

    .. code-block:: matlab

       cond(11).name = 'Awesome tissue';
       cond(11).value = 10;

    \
