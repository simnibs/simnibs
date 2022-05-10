.. _faq:

FAQ
===



Can SimNIBS simulate TACS?
--------------------------

Yes! Frequencies used in TACS are usually low, so we can employ a *quasi-static approximation*. This means that if we have a sinusoidal current with a frequency :math:`f` and amplitude :math:`I_0`

.. math::

   I(t) = I_0\sin\left(2 \pi f t\right)

The electric field will vary in time and with the same frequency of the input current and with the same phase:

.. math::

   \boldsymbol E (t) = \boldsymbol E_0\sin\left(2 \pi f t\right)

Where :math:`\boldsymbol E_0` is the electric field obtained with the current :math:`I_0`. In practice, this means you can just simulate the electric field obtained at peak currents, and the temporal variations of the current will just scale the field, and not change its distribution in the brain.

When multiple channels are involved, and especially when they are out of phase, more care needs to be taken. Please see `(Saturnino et. al., 2017) <https://doi.org/10.1016/j.neuroimage.2017.09.024>`_ for a more detailed discussion



What is the dI/dt in TMS simulations?
-------------------------------------
:math:`dI/dt` is the speed of variation of the current through the coil. This value depends on the coil model, the stimulator model and the pulse intensity. In some stimulators, it is displayed on the screen after a pulse. As :math:`dI/dt` varies in time, usually the value at the beginning of the pulse is taken, which corresponds to the peak :math:`dI/dt` for most pulse shapes.

:math:`dI/dt` depends approx. linearly on the pulse intensity set by the user, meaning that a setting of 80% maximum stimulator output (MSO) will also give 80% of the maximal :math:`dI/dt`. That is, :math:`^{dI}/_{dt}=0.8*{^{dI}/_{dt}}_{max}` when stimulating at 80% MSO.

For the 25 coil models in the subfolder Drakaki_BrainStim_2022, the :math:`{^{dI}/_{dt}}_{max}` values are listed in Table 2 of the the `corresponding paper <https://doi.org/10.1016/j.brs.2022.04.017>`_ for the most commonly used biphasic stimulators. They can be used together with the %MSO to calculate the :math:`dI/dt` for the simulations.

The Electric Field varies linearly with :math:`^{dI}/_{dt}`


.. math::

  \boldsymbol E(t) = \boldsymbol E_0 \frac{dI}{dt}(t)


Where :math:`\boldsymbol E_0` is the electric field obtained with a unit :math:`^{dI}/_{dt}` value.




Are the SimNIBS coordinates in MNI Space?
------------------------------------------

No. SimNIBS uses coordinates defined in the subject space. But we provide many tools to transform between subject and MNI spaces. Please see :ref:`coords_doc` for more information.



Electric Field Magnitude? Normal?
---------------------------------

The Electric Field is a `vector field <https://www.khanacademy.org/math/multivariable-calculus/thinking-about-multivariable-function/ways-to-represent-multivariable-functions/a/vector-fields>`_, meaning it has a direction everywhere in space.

Vectors are hard to visualize (even though you can do it in SimNIBS!), so we often use scalars to represent information about the vectors. The one we use the most is the **magnitude** (short: *magn*), which is the *length* or *strength* of the vector, irrespective of direction (Note: up to version 3, the *magnitude* was called *norm*, as mathematicllay it is the vector norm).

Another useful quantity is the **normal** with respect to a given surface, normally the middle cortical surface. The **normal** gives us the field that is incoming or outgoing from the surface. 

.. figure:: images/tutorial_normal.png
   :align: center
   
   From `Antonenko et. al. (2019) <https://doi.org/10.1016/j.brs.2019.03.072>`_

\

The effect that the different polarities have on neuronal modulation suggests that the field **normal** is important in tDCS. For TMS, it is harder to determine which field component if any is more relevant for stimulation.


Which Units does SimNIBS Use?
-------------------------------

SimNIBS almost always uses units in the International System of Units (SI units).
The exceptions are head model node positions, which are in millimiter and current values in the GUI, which are in mA.

This means we have the units:

.. list-table::
   :widths: 30 10
   :header-rows: 1

   * - Quantity
     - Units
   * - Electric field (and respective magnitude, normal, ...)
     - V/m
   * - Current Density (and respective magnitude, normal, ...)
     - A/m²
   * - Conductivities
     - S/m
   * - Electrode Currents
     - mA (GUI) / A (scripts)
   * - Coil dI/dt
     - A/s
   * - Mesh node positions
     - mm
   * - Magnetic Field
     - T

