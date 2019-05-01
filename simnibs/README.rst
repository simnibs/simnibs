SimNIBS
#######

Simulation of Non-Invasive Brain Stimulation

Python 2.7 library for simulations of TMS and tDCS electric fields www.simnibs.org

.. code-block::python
	>>> import simnibs
	# Run TMS simulation
	>>> tms_sim = simnibs.tms_simulation('almi5.msh', 'Magstim_70mm_Fig8', 'AP',
	>>>                                  [-42.1, -13.8, 54.3], 'tms_simulation')
	>>> tms_sim.run()
	>>> tms_sim.show_result()
	# Run tDCS simulation   	
        >>> stimulation_el = simnibs.tdcs_electrode('rectangular', [5, 5], [-42.1, -13.8, 54.3], 1.0)
        >>> return_el = simnibs.tdcs_electrode('rectangular', [5, 7], [19.1, 80.0, 2.4], -1.0)
	>>> tdcs_sim = simnibs.tdcs_simulation('almi5.msh', [stimulation_el, return_el], 'tdcs_simulation')

