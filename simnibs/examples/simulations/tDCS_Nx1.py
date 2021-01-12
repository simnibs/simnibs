"""
Example on running a SimNIBS tDCS simulation for a Nx1 montage in Python.
    
    Run with:
        simnibs_python tDCS_Nx1.py

    Copyright (C) 2020 Axel Thielscher
"""
from copy import deepcopy
import numpy as np
from scipy.spatial.transform import Rotation as R
import os
from simnibs import sim_struct, run_simnibs, mesh_io, file_finder


def sphereFit(pts, bounds = None):
    """
    Fit a circle or sphere to a point cloud.
    
    returns the radius and center points of the best fit sphere
    (adapted from https://jekel.me/2015/Least-Squares-Sphere-Fit/)
    
    Parameters
    ----------
    pts : array (Nx2) or (Nx3)
        point cloud
        
    Returns
    -------
    R : float64
        radius
    centre: ndarray (2,) or (3,)
        centre position
        
    """
    A = np.hstack((2*pts, np.ones((pts.shape[0],1)) ))
    f = np.sum(pts**2,1)
    C, residuals, rank, singval = np.linalg.lstsq(A,f,rcond=None)
 
    dim = pts.shape[1] 
    R = np.sqrt(np.sum(C[0:dim]**2) + C[dim])
    return R, C[0:dim]


def sph2cart(az, el, r): # phi, theta, radius
    """Conversion from spherical to cartesian coordinates."""
    rcos_theta = r * np.cos(el)
    pts = np.zeros(( (3,) + rcos_theta.shape ))
    pts[0,:] = rcos_theta * np.cos(az)
    pts[1,:] = rcos_theta * np.sin(az)
    pts[2,:] = r * np.sin(el)
    return pts


def show_for_debugging(m,sph_centre,radius,P_centre,P_surround,surround_fit,M_sph):
    """Show some results in gmsh for debugging."""
    import tempfile
    fn_geo=tempfile.NamedTemporaryFile(suffix='.geo').name    
    mesh_io.write_geo_spheres(sph_centre.reshape((1,3)),
                              fn_geo, name = 'center', mode = 'bw')
    mesh_io.write_geo_spheres(( M_sph @ [radius,0,0,1] )[:3].reshape((1,3)),
                              fn_geo, name = 'x', mode = 'ba')
    mesh_io.write_geo_spheres(( M_sph @ [0,radius,0,1] )[:3].reshape((1,3)),
                              fn_geo, name = 'y', mode = 'ba')
    mesh_io.write_geo_spheres(( M_sph @ [0,0,radius,1] )[:3].reshape((1,3)),
                              fn_geo, name = 'z', mode = 'ba')
    
    TH_DEBUG = np.arange(-1.0,1.01,0.1)*np.pi
    PHI_DEBUG = np.arange(0.,1.01,0.05)*2*np.pi
    TH_DEBUG, PHI_DEBUG = np.meshgrid(TH_DEBUG, PHI_DEBUG)
    TH_DEBUG = TH_DEBUG.flatten()
    PHI_DEBUG = PHI_DEBUG.flatten()
    R_DEBUG = radius*np.ones_like(TH_DEBUG)
    
    pts=sph2cart(PHI_DEBUG, TH_DEBUG, R_DEBUG)
    pts = np.vstack(( pts, np.ones((1,pts.shape[1])) ))
    mesh_io.write_geo_spheres((M_sph @ pts)[:3,:].T, fn_geo,
                              name = 'sphere', mode = 'ba')
    
    mesh_io.write_geo_spheres( P_centre.reshape((1,3)),
                               fn_geo, name = 'centre', mode = 'ba')
    for i in range(len(P_surround)):
        mesh_io.write_geo_spheres( P_surround[i].reshape((1,3)),
                                  fn_geo, name = 'surr '+str(i), mode = 'ba')
    
    N_pts = 50
    for i in range(len(surround_fit)):
        tmp_centre = surround_fit[i][0]
        tmp_r = surround_fit[i][1]
        tmp_theta = surround_fit[i][2]
        tmp_theta_z0 = surround_fit[i][3]
        tmp_M = surround_fit[i][4]
        
        tmp_arc = np.vstack(( 
            tmp_r*np.sin(tmp_theta_z0 + (tmp_theta-tmp_theta_z0)*np.arange(N_pts)/(N_pts-1)) + tmp_centre[0],
            np.zeros((1,N_pts)),
            tmp_r*np.cos(tmp_theta_z0 + (tmp_theta-tmp_theta_z0)*np.arange(N_pts)/(N_pts-1)) + tmp_centre[1],
            np.ones((1,N_pts))
            ))
        tmp_arc=(tmp_M @ tmp_arc)[:3].T
        mesh_io.write_geo_spheres(tmp_arc, fn_geo, name = 'arc '+str(i), mode = 'ba')
          
    vis = mesh_io.gmsh_view.Visualization(m)
    vis.add_merge(fn_geo)
    vis.show()
    os.remove(fn_geo)


def get_surround_pos(center_pos, fnamehead, radius_surround = 50, N = 4, 
                     pos_dir_1stsurround = None):
    """
    Determine the positions of surround electrode.

    Parameters
    ----------
    center_pos : array (3,) or string
        Center position of the central electrode.
    fnamehead : string
        Filename of head mesh.
    radius_surround : float, optional
        Distance (centre-to-centre) between the centre and 
        surround electrodes (mm). The default is 50.
    N : int, optional
        Number of surround electrodes. The default is 4.
    pos_dir_1stsurround : array (3,) or string, optional
        A position indicating the direction from center_pos to 
        the position of the first surround electrode. The default is None.

    Returns
    -------
    P_surround : list of arrays (3,)
        List of the centre positions of the surround electrodes.

    """
    DEBUG=False
    if DEBUG:
        import matplotlib.pyplot as plt
        from matplotlib import cm
        cmap = cm.get_cmap('hsv', N)
        cmap = cmap(np.linspace(0, 1, N)) 
        
    # replace electrode name with position if needed
    # and get skin ROI around centre position
    ff = file_finder.SubjectFiles(fnamehead = fnamehead)
    tmp = sim_struct.ELECTRODE()
    tmp.centre = center_pos
    tmp.substitute_positions_from_cap(ff.get_eeg_cap())
    
    m = mesh_io.read_msh(fnamehead)
    idx = (m.elm.elm_type == 2)&( (m.elm.tag1 == 1005) | (m.elm.tag1 == 5) )
    m = m.crop_mesh(elements = m.elm.elm_number[idx])
    P_centre = m.find_closest_element(tmp.centre)
    idx = np.sum((m.nodes[:] - P_centre)**2,1) <= (radius_surround+10)**2
    m = m.crop_mesh(nodes = m.nodes.node_number[idx])
    
    
    # fit sphere to skin ROI to build local coordinate system
    #   origin: sphere center
    #   x-axis: direction of first surround
    #   z-axis: from sphere center to centre electrode
    r_sph, sph_centre = sphereFit(m.nodes[:])
    theta_on_sph = radius_surround/r_sph
    
    M_sph = np.eye(4)
    M_sph[:3,3] = sph_centre
    tmp = P_centre - sph_centre
    M_sph[:3,2] = tmp/np.linalg.norm(tmp)
    # direction of first surround
    if pos_dir_1stsurround is not None:
        # replace electrode name with position if needed
        tmp = sim_struct.ELECTRODE()
        tmp.centre = pos_dir_1stsurround
        tmp.substitute_positions_from_cap(ff.get_eeg_cap())
        tmp = tmp.centre - P_centre # this is not orthogonal to Z
    else:
        # get a vector orthogonal to z-axis   
        tmp = np.cross(M_sph[:3,2],np.eye(3))
        tmp = tmp[:, np.argmax( np.linalg.norm(tmp,axis=1) )]
    M_sph[:3,1] = np.cross(M_sph[:3,2], tmp)
    M_sph[:3,1] /= np.linalg.norm(M_sph[:3,1])
    M_sph[:3,0] = np.cross(M_sph[:3,1], M_sph[:3,2])
    
    
    # fit arcs to the skin to get the distances accurate
    N_pts = 50
    arc = np.vstack(( r_sph*np.sin(theta_on_sph*np.arange(N_pts)/(N_pts-1)),
                      np.zeros((1,N_pts)),
                      r_sph*np.cos(theta_on_sph*np.arange(N_pts)/(N_pts-1)),
                      np.ones((1,N_pts)) ))
    P_surround = []
    surround_fit = []
    for phi in range(N):
        M_rot = np.eye(4)
        M_rot[:3,:3] = R.from_euler('z', phi/N*2*np.pi ).as_dcm()
        M_to_world = M_sph @ M_rot
        M_from_world = np.linalg.inv(M_to_world)
        
        # project skin points into XZ-plane that contains the arc
        Pts = (M_to_world @ arc).T
        Pts[:,:3] = m.find_closest_element(Pts[:,:3])
        Pts = M_from_world @ Pts.T
        
        # fit individual arc
        r_arc, arc_centre = sphereFit(Pts[(0,2),:].T)
        
        if np.abs(arc_centre[0]) > r_arc:
            # z-axis does not intersect with circle 
            # --> use initial sphere instead            
            r_arc = r_sph
            arc_centre *= 0
            
        theta_z0_on_arc = -np.arcsin(arc_centre[0]/r_arc)
        if arc_centre[1]<np.mean(Pts[2,:]):
            theta_on_arc = radius_surround/r_arc + theta_z0_on_arc
        else:
            # best fitting arc has opposite curvature compared
            # to initial sphere
            theta_z0_on_arc = - theta_z0_on_arc + np.pi
            theta_on_arc = theta_z0_on_arc - radius_surround/r_arc
            
        # get centre of surround electrode
        tmp = np.array(( r_arc*np.sin(theta_on_arc) + arc_centre[0],
                         0.,
                         r_arc*np.cos(theta_on_arc) + arc_centre[1],
                         1. ))    
        P_surround.append(m.find_closest_element( (M_to_world @ tmp).T[:3] ))
        surround_fit.append((arc_centre, r_arc, theta_on_arc, theta_z0_on_arc, M_to_world))
    
        if DEBUG:
            arc_fit = np.vstack(( 
                 r_arc*np.sin(theta_z0_on_arc + (theta_on_arc-theta_z0_on_arc)*np.arange(N_pts)/(N_pts-1)) + arc_centre[0],
                 r_arc*np.cos(theta_z0_on_arc + (theta_on_arc-theta_z0_on_arc)*np.arange(N_pts)/(N_pts-1)) + arc_centre[1]
                 ))
            plt.plot(Pts[0,:],Pts[2,:],marker='o',color=cmap[phi,:],linestyle='None')
            plt.plot(arc_fit[0,:],arc_fit[1,:],color=cmap[phi,:])
    if DEBUG:
        plt.plot(arc[0,:],arc[2,:],linewidth=4,color='black')    
        plt.show()
        print(np.linalg.norm(np.array(P_surround)-P_centre,axis=1))
        show_for_debugging(m,sph_centre,r_sph,P_centre,P_surround,surround_fit,M_sph)

    return P_surround


def expand_to_center_surround(S, fnamehead, radius_surround = 50, N = 4, 
                              pos_dir_1stsurround = None, multichannel = False):
    """
    Generate a center-surround montage (standard: 4x1) from a TDCSLIST.
    
    The TDCSLIST has to contain only the center electrode. Copies of this
    electrode are then placed in a circle around the centre 

    Parameters
    ----------
    S : TDCSLIST
        TDCSLIST with the center electrode.
    fnamehead : string
        Filename of the head mesh.
    radius_surround : float, optional
        Distance (centre-to-centre) between the centre and surround
        electrodes. The default is 50.
    N : int, optional
        Number of surround electrodes. The default is 4.
    pos_dir_1stsurround : array (3,) or string, optional
        A position indicating the direction from center_pos to 
        the position of the first surround electrode. The default is None.
    multichannel : Boolean, optional
        When set to True, a multichannel stimulator with each suround channel 
        receiving 1/N-th of the of the center channel will be simulated
        (standard: False, i.e. all surround electrodes connected to the 
         same return channel).

    Returns
    -------
    S : TDCSLIST
        TDCSLIST with the surround electrodes added.

    """
    if S.type != 'TDCSLIST':
        raise TypeError('The first parameter needs to be a TDCSLIST.')
    if len(S.electrode) != 1:
        raise ValueError('The TDCSLIST has to contain exactly one ELECTRODE.')
    if not os.path.isfile(fnamehead):
        raise IOError('Could not find head mesh: {0}'.format(fnamehead))
    
    C = S.electrode[0]
    C.channelnr = 1  # Connect center to channel 1
    if not len(C.name):
        C.name = 'centre'
        
    # set surround channels and current strengths
    if type(S.currents) == float:
        C_current = S.currents
    else:
        C_current = S.currents[0]
        
    if multichannel:
        S.currents = -C_current/N*np.ones(N+1)
        S.currents[0] = C_current
        Channel_surround = np.arange(2,N+2)
    else:
        S.currents = [C_current, -C_current]
        Channel_surround = 2*np.ones(N,dtype = int)
    
    # get centers of surround electrodes
    P_surround = get_surround_pos(C.centre, fnamehead, radius_surround = radius_surround,
                                  N = N, pos_dir_1stsurround = pos_dir_1stsurround)
    
    # get direction vector
    ydir = []
    if len(C.pos_ydir):    
        ff = file_finder.SubjectFiles(fnamehead = fnamehead)
        tmp = deepcopy(C)
        tmp.substitute_positions_from_cap(ff.get_eeg_cap())
        ydir = tmp.pos_ydir - tmp.centre
    
    # add surround electrodes to TDCSLIST
    for i in range(N):
        S.electrode.append(deepcopy(C))
        El = S.electrode[-1]
        El.centre = P_surround[i]
        El.channelnr = Channel_surround[i]
        El.name = 'surround '+str(i+1)
        if len(ydir):
            El.pos_ydir = El.centre + ydir
    return S


###     SETUP
###################
S = sim_struct.SESSION()
S.fnamehead = 'ernie.msh'  # head mesh
S.pathfem = 'tdcs_Nx1' # output directory for the simulation
#S.map_to_surf = True # map to subject's middle gray matter surface (optional)

tdcs_list = S.add_tdcslist()
tdcs_list.currents = 0.001  # Current flow through center channel (mA)

# define the center electrode
center = tdcs_list.add_electrode()
center.centre = 'C3'  # Place it over C3
center.shape = 'ellipse'  # round shape
center.dimensions = [30, 30]  # 30 mm diameter
center.thickness = [2, 1]  # 2 mm rubber electrodes on top of 1 mm gel layer

# parameters for setting up the surround electrodes
radius_surround = 80 # distance (centre-to-centre) between the centre and 
                     # surround electrodes (optional; standard: 50 mm)
pos_dir_1stsurround = 'C4' # a position indicating the direction in which the 
                           # first surround electrode should be placed 
                           # (optional; standard: None)
N = 4 # number of surround electrodes (optional; standard: 4)
multichannel = True # when set to True: Simulation of multichannel stimulator 
                     # with each suround channel receiving 1/N-th of the
                     # center channel (optional; standard: False, i.e. all 
                     # surround electrodes connected to the same channel)

# set up surround electrodes
tdcs_list = expand_to_center_surround(tdcs_list, S.fnamehead, radius_surround, 
                                      N, pos_dir_1stsurround, multichannel)
    

### RUN SIMULATION
###################
run_simnibs(S)
