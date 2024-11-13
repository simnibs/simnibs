"""
Example for calculating TI envelopes from pre-calculated leadfields

Copyright (c) 2022 SimNIBS developers. Licensed under the GPL v3.
"""

import copy
import numpy as np

from simnibs import mesh_io
from simnibs.utils import TI_utils as TI

# load lead field       
leadfield_hdf = 'leadfield/ernie_leadfield_EEG10-10_UI_Jurak_2007.hdf5'
leadfield, mesh, idx_lf = TI.load_leadfield(leadfield_hdf) 


# specify two electrode pairs with current intensities
TIpair1 = ['F5','P5',0.001]
TIpair2 = ['F6','P6',0.001]

# get fields for the two pairs
ef1 = TI.get_field(TIpair1,leadfield,idx_lf)
ef2 = TI.get_field(TIpair2,leadfield,idx_lf)

# add to mesh for later visualization
mout = copy.deepcopy(mesh)
hlpvar = mesh_io.NodeData(ef1,mesh=mout)
mout.add_node_field(hlpvar.norm(),'E_magn1')
hlpvar = mesh_io.NodeData(ef2,mesh=mout)
mout.add_node_field(hlpvar.norm(),'E_magn2')


# Option 1: get maximal TI amplitude
TImax = TI.get_maxTI(ef1,ef2)
mout.add_node_field(TImax,'TImax') # for visualization


# Option 2: get TI amplitudes along x, y and z
TIamp = TI.get_dirTI(ef1,ef2,[1,0,0])
mout.add_node_field(TIamp,'TIamp_x') # for visualization

TIamp = TI.get_dirTI(ef1,ef2,[0,1,0])
mout.add_node_field(TIamp,'TIamp_y') # for visualization

TIamp = TI.get_dirTI(ef1,ef2,[0,0,1])
mout.add_node_field(TIamp,'TIamp_z') # for visualization


# Option 3: get TI amplitudes along local normal orientation
surf_normals = mesh.nodes_normals().value
TIamp = TI.get_dirTI(ef1,ef2,surf_normals)
mout.add_node_field(TIamp,'TIamp_localnorm') # for visualization


mesh_io.write_msh(mout,'TI_via_leadfields.msh')
v = mout.view(
    visible_tags=[1, 2, 1006],
    visible_fields='TImax',    
    )
v.write_opt('TI_via_leadfields.msh')
mesh_io.open_in_gmsh('TI_via_leadfields.msh', True)

