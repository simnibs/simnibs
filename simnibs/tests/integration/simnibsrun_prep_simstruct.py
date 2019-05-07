import argparse
import sys
import os
import simnibs.simulation.sim_struct as sim_struct
import numpy as np

try:
    path_simnibs = os.path.realpath(os.environ["SIMNIBSDIR"])
except:
    raise IOError("SIMNIBSDIR environment variable not found. It should be set inside ~/.bashrc")


def parse_arguments(argv):
    parser = argparse.ArgumentParser(prog="testing_prep_simstruct.py",
                                     description="Reads a mesh file, a CSV file with"
                                     " coil and electrode positions, and a prototype"
                                     " simstruct; creates a ready-to-be-used simstruct")
    parser.add_argument("subid",
                        help="Subject ID")
    parser.add_argument('position_csv',
                        help="CSV file with electrode and coil positions")
    parser.add_argument('prototype_structure',
                        help="Prototype structure with simulation information")
    return parser.parse_args(argv)


if __name__ == '__main__':
    args = parse_arguments(sys.argv[1:])

    coils = (os.path.join(path_simnibs,'ccd-files','MagVenture_MC_B70.ccd') ,
             os.path.join(path_simnibs,'ccd-files','MagVenture_MC_B70.nii.gz'))
    
    S = sim_struct.read_matlab_sim_struct(args.prototype_structure)
    S.subpath = 'm2m_' + args.subid
    S.fields='eEjJ'
    S.pathfem = 'simu'
    S.fname_tensor = sim_struct.get_dir_tensor_from_m2m_folder(S.subpath)
    S.fnamehead = sim_struct.get_mesh_name_from_m2m_folder(S.subpath)
    S.fiducials.from_csv(args.position_csv)
    
    if (os.path.isfile(os.path.join(S.subpath,'mri2mesh_log.html')) or
        os.path.isfile(os.path.join(S.subpath,'mask_prep','MASK_CAT_GM.nii.gz'))):
        S.map_to_surf=True
        S.map_to_fsavg=True
        S.map_to_vol=True
        S.map_to_MNI=True

    for i in range(0,len(S.poslists)):
        if S.poslists[i].type.upper() == 'TDCSLIST':
            print str(i) + ' TDCSLIST'
            S.poslists[i].add_electrodes_from_csv(args.position_csv)
        elif S.poslists[i].type.upper() == 'TMSLIST':
            print str(i) + ' TMSLIST'
            S.poslists[i].pos=[]
            S.poslists[i].add_positions_from_csv(args.position_csv)
            if i >= len(coils):
                S.poslists[i].fnamecoil=coils[0]
            else:
                S.poslists[i].fnamecoil=coils[i]
    
    S.prepare()
    sim_struct.save_matlab_sim_struct(S, args.subid + '_test_simstruct.mat')
