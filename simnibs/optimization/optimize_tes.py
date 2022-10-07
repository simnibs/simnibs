import copy
import csv
import re
import os
import time
import glob
import functools
import logging
import gc

import numpy as np
import scipy.spatial
import h5py
import nibabel

from . import optimize_tms
from . import optimization_methods
from . import ADMlib
from ..simulation import fem
from ..simulation import cond
from ..simulation.sim_struct import SESSION, TMSLIST, SimuList, save_matlab_sim_struct
from ..mesh_tools import mesh_io, gmsh_view
from ..utils import transformations
from ..utils.simnibs_logger import logger
from ..utils.file_finder import SubjectFiles
from ..utils.matlab_read import try_to_read_matlab_field, remove_None


# TODO: implement new tes fast optimize class
# TODO: adapt imports
# TODO: adapt matlab_tools/opt_struct.m and include new class

# TODO: Kodierung der Kanäle beim Array Layout! Welche Kanäle zusammen geschaltet sind
# 1:1 Mapping der Elektroden -> von Neumann Kanal 1 -> Sink 1, Kanal 2 -> Sink 2
# 1 Kanal -> mehrere sink elektroden -> fake Dirichlet -> Kanal 1 -> Sink 1,2,3 (zusammengeschaltet) gleiche Spannung


class TDCSoptimize():
    ''' Defines a fast tdcs optimization problem

    Parameters
    --------------
    leadfield_hdf: str (optional)
        Name of file with leadfield
    max_total_current: float (optional)
        Maximum current across all electrodes (in Amperes). Default: 2e-3
    max_individual_current: float (optional)
        Maximum current for any single electrode (in Amperes). Default: 1e-3
    max_active_electrodes: int (optional)
        Maximum number of active electrodes. Default: no maximum
    name: str (optional)
        Name of optimization problem. Default: optimization
    target: list of TDCStarget objects (optional)
        Targets for the optimization. Default: no target
    avoid: list of TDCSavoid objects
        list of TDCSavoid objects defining regions to avoid


    Attributes
    --------------
    leadfield_hdf: str
        Name of file with leadfield
    max_total_current: float (optional)
        Maximum current across all electrodes (in Amperes). Default: 2e-3
    max_individual_current: float
        Maximum current for any single electrode (in Amperes). Default: 1e-3
    max_active_electrodes: int
        Maximum number of active electrodes. Default: no maximum

    ledfield_path: str
        Path to the leadfield in the hdf5 file. Default: '/mesh_leadfield/leadfields/tdcs_leadfield'
    mesh_path: str
        Path to the mesh in the hdf5 file. Default: '/mesh_leadfield/'

    The two above are used to define:

    mesh: simnibs.msh.mesh_io.Msh
        Mesh with problem geometry

    leadfield: np.ndarray
        Leadfield matrix (N_elec -1 x M x 3) where M is either the number of nodes or the
        number of elements in the mesh. We assume that there is a reference electrode

    Alternatively, you can set the three attributes above and not leadfield_path,
    mesh_path and leadfield_hdf

    lf_type: None, 'node' or 'element'
        Type of leadfield.

    name: str
        Name for the optimization problem. Defaults tp 'optimization'

    target: list of TDCStarget objects
        list of TDCStarget objects defining the targets of the optimization

    avoid: list of TDCSavoid objects (optional)
        list of TDCSavoid objects defining regions to avoid

    open_in_gmsh: bool (optional)
        Whether to open the result in Gmsh after the calculations. Default: False

    Warning
    -----------
    Changing leadfield_hdf, leadfield_path and mesh_path after constructing the class
    can cause unexpected behaviour
    '''

    def __init__(self, leadfield_hdf=None,
                 max_total_current=2e-3,
                 max_individual_current=1e-3,
                 max_active_electrodes=None,
                 name='optimization/tdcs',
                 target=None,
                 avoid=None,
                 open_in_gmsh=True,
                 electrode=None):
        self.leadfield_hdf = leadfield_hdf
        self.max_total_current = max_total_current
        self.max_individual_current = max_individual_current
        self.max_active_electrodes = max_active_electrodes
        self.leadfield_path = '/mesh_leadfield/leadfields/tdcs_leadfield'
        self.mesh_path = '/mesh_leadfield/'
        self.open_in_gmsh = open_in_gmsh
        self._mesh = None
        self._leadfield = None
        self._field_name = None
        self._field_units = None
        self.name = name
        # I can't put [] in the arguments for weird reasons (it gets the previous value)
        if target is None:
            self.target = []
        else:
            self.target = target
        if avoid is None:
            self.avoid = []
        else:
            self.avoid = avoid

        # TODO: precompute all node areas
        self.areas = self.get_node_areas()

        # TODO: create array_layout class with electrodes
        self.array_layout = array_layout

        # TODO: write postproc prepare function calls the new stuff from fem.py see SURFACE but rename to RegionOfInterest, creates sF matrix for SPR etc...
        self.postproc_prepare()

        # determine point indices where the electrodes may be applied during optimization
        self.valid_skin_region()

        # gauge problem (find node in COG of head and move to end), modifies self.mesh
        self.gauge_mesh()

        # TODO: solver_options: use PARDISO solver as standard
        # assemble FEM matrix
        self.fem=FEMSystem.tdcs_neumann(
            self.mesh, elm_cond, self.mesh.node.nr,  # last element is 0 was shifted from center
            solver_options=solver_options,
            input_type='node'
        )

    def get_node_areas(self):
        """
        Computes associated areas of nodes
        :return:
        """
        self.areas = 1

    def gauge_mesh(self):
        """

        :return:
        """
        pass

    def postproc_prepare(self):
        """
        Prepares postprocessing by computing sF matrices for each ROI
        :return:
        """
        pass

    def valid_skin_region(self):
        """
        Computes the node indices on the scalp surface where the electrode can be applied
        """
        pass

    def update_electrode(self, location_parameters):
        """
        Updates
        :return:
        """
        # this function is currently in sim_struct and rather slow
        self.array_layout.get_surround_pos(center_pos, fnamehead, radius_surround=50, N=4,
                         pos_dir_1stsurround=None, phis_surround=None,
                         tissue_idx=1005, DEBUG=False)
        self.array_layout_node_idx = 1

    def update_rhs(self):
        """
        Update RHS with new electrode positions
        :return:
        """
        self.rhs = self.fem.assemble_tdcs_neumann_rhs([np.hstack((self.array_layout_node_idx))], [np.hstack((I))], input_type='node', areas=self.areas)


    # TODO: From Fang (TMS) adapt accordingly
    # TODO: surface can be a ROI class with a .get_fields(v) method, is initialized with either "TMS" (with dAdt or "TES" raw)
    # TODO: change dataType ROI specific in ROI class init
    def update_field(self, matsimnibs, coil_id, amplitude, surfaceIDs, dataType):
        '''Calculate the E field with the different coil positions. The input matsimnibs is a 4x4 affine matrix which defines the coil position in 3D'''

        start = time()
        a_affine = self.coil[coil_id].a_affine
        a_field = self.coil[coil_id].a_field
        didtmax = self.coil[coil_id].didtmax

        assert dataType in [0, 1]



        start = time()

        # update electrode positions
        self.update_electrode()

        # Assemble RHS
        self.update_rhs()

        logger.debug('Assemble the right hand side force vector: {0:.6f}s'.format(time() - start))

        start = time()

        # Solve the linear equation Mx = F. M is the stiffness matrix and F is the force vector.
        x = self.solve_and_normalize(rhs)

        # TODO: generalize
        # Dirichlet node in case of TMS (last one set to zero)
        v = np.append(x, 0)

        logger.debug('Estimate the x: {0:.2f}s'.format(
            time() - start))

        start = time()

        m = []

        # TODO: this is done in the new ROI class in the .calc_fields() method, TMS needs dAdt !!!
        # loop over the surfaces
        for s in surfaceIDs:

            # Post-processing: estimate the E field
            # step1: get the E field in all tetrahedra
            # v: (number_of_nodes, 1)
            fields = np.einsum('ijk,ij->ik', self.surface[s].gradient, - (v * 1e3)[self.surface[s].node_index_list])

            if dataType == 0:
                # calculate the magnitutde of E
                fields = np.linalg.norm(fields, axis=1, keepdims=True)

            m.append(self.surface[s].sF @ fields)

        logger.debug('Estimate the E field on the surfaces: {0:.2f}s'.format(
            time() - start))

        return m

    def solve_and_normalize(self, node_idx, I, v_norm, I_norm):
        # set RHS (in fem.py, check for speed)
        b = self.fem.assemble_tdcs_neumann_rhs([np.hstack((self.array_layout_node_idx))], [np.hstack((I))], input_type='node')
        # solve
        v = self.fem.solve(b)
        v_elec = [v[self.array_layout_node_idx[k] - 1] for k in range(len(I))]

        # TODO: add postprocessing step from Fang to get e-fields

        v_mean = [np.mean(v_elec[k]) for k in range(len(I))]
        for k in range(len(I)):
            vn = (v_elec[k] - np.mean(v_elec[k])) / np.std(v_mean)
            In = (I[k] - np.mean(I[k])) / np.mean(I[k])

            v_norm[k] = np.append(v_norm[k], vn.reshape(1, len(vn)), axis=0)
            I_norm[k] = np.append(I_norm[k], In.reshape(1, len(In)), axis=0)
        return v, v_elec, v_norm, I_norm

    def run(self, cpus=1, allow_multiple_runs=False, save_mat=True, return_n_max=1):
        ''' Runs the optimization problem

        Parameters
        -------------
        fn_out_mesh: str
            If set, will write out the electric field and currents to the mesh

        fn_out_mesh: str
            If set, will write out the currents and electrode names to a CSV file


        Returns
        ------------
        currents: N_elec x 1 ndarray
            Optimized currents. The first value is the current in the reference electrode
        '''
        pass
