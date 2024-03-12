import os
import pytest
import simnibs
import numpy as np
from ... mesh_tools import mesh_io, surface
from ... import SIMNIBSDIR
from .. opt_struct import TESoptimize, valid_skin_region
from ...utils.file_finder import Templates


class TestRegionOfInterest:
    # Not possible without m2m_* folder because of mapping from MNI to subject space
    # def test_valid_skin_region(self):
    #     fn_mesh = os.path.join(SIMNIBSDIR, '_internal_resources', 'testing_files', 'sphere3.msh')
    #     mesh = mesh_io.read_msh(fn_mesh)
    #     skin_surface = surface.Surface(mesh=mesh, labels=1005)
    #     fn_electrode_mask = Templates().mni_volume_upper_head_mask
    #     skin_surface_valid = valid_skin_region(skin_surface, mesh, fn_electrode_mask, additional_distance=0)

    # Not possible without m2m_* folder
    # def test_relabel_internal_air(self):
    #     fn_mesh = os.path.join(SIMNIBSDIR, '_internal_resources', 'testing_files', 'sphere3.msh')
    #     mesh = mesh_io.read_msh(fn_mesh)
    #     mesh_relabel = relabel_internal_air(m=mesh,
    #                                         subpath=os.path.split(self.mesh.fn)[0],
    #                                         label_skin=1005,
    #                                         label_new=1099,
    #                                         label_internal_air=501)

    def test(self):
        pass
