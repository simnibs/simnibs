function S = opt_struct(type)
%
% create an empty data structure
% to set up a simnibs optimization
%
% S = opt_struct(type)
% 
% Guilherme Saturnino, 2019

%    This program is part of the SimNIBS package.
%    Please check on www.simnibs.org how to cite our work in publications.
%
%    Copyright (C) 2020 Guilherme Saturnino
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.


validtypes={'TMSoptimize', 'TDCSoptimize', 'TDCStarget', 'TDCSavoid', 'TDCSDistributedOptimize', ...
            'TmsFlexOptimization', 'RegionOfInterest', 'TesFlexOptimization', 'CircularArray', ...
            'ElectrodeArrayPair'};

if ~any(strcmp(type,validtypes))
    disp(validtypes)
    error('structure type has to be one of the above')
end

S.type=type;

switch S.type
    case 'TMSoptimize'
       S.fnamehead=''; % same as ${subID}.msh created by mri2mesh or headreco
       S.subpath = ''; % path to the 'm2m_{subID}' folder created by mri2mesh or headreco (OPTIONAL, filled from fnamehead)
       S.pathfem = 'tms_optimization/';   % path to save the results (OPTIONAL)
       S.fname_tensor = ''; % file name of the diffusion tensors (OPTIONAL, filled from fnamehead)
       S.fnamecoil =  '';      % to chose from inside coil_models
       S.cond = standard_cond;   % list of conductivities
       S.anisotropy_type = 'scalar'; % can be 'scalar' (use isotropic values), 'dir' (direct mapping),'mc' (mean conductivity from direct mapping),'vn' (volume normalized); optional
       S.aniso_maxratio = 10; % maximal ratio between largest eigenvalue and the two other eigenvalues of conductivity tensor
       S.aniso_maxcond = 2; % maximal directional conductivity in [S/m] (i.e. max eigenvalue of conductivity tensor)
       S.tissues = [2]; % list, tissues where the target is defined (Optional)
       S.target = []; % Position of the optimization target, in head coordinates
       S.target_direction = []; % Direction of electric field to be optimized (Optional).
       S.target_size = 5; % Size of target, in mm
       S.centre = []; % Position in scalp to use as a reference for the search space (Optional) . By default, will project the target to the scalp.
       S.pos_ydir = []; % Reference position for the coil Y axis, with respect to the target (or the pos variable, if it is defined). If left empty, will search positions in a 360 degrees radius.
       S.distance = 4; % Distance between coil and scalp (in mm)
       S.didt = 1e6; % Coil current (in A/m)
       S.search_radius = 20; % Radius around the "pos" to search, in mm
       S.spatial_resolution = 5; % Spatial resolution for serach, in mm
       S.search_angle = 360; % Range of angles to search (in degrees)
       S.angle_resolution = 30; % Resolution to use for angles (in degrees)
       S.open_in_gmsh = true; % Wether to open simulation result in Gmsh
       S.solver_options = ''; % FEM solver options
       S.method = 'direct'; % Solution method, either 'direct' or 'ADM'. The former is only valid with .ccd coil files or .tcd files that only contain dipole elements

    case 'TDCSoptimize'
        S.leadfield_hdf = ''; % Name of HDF5 file with leadfield
        S.name = 'optimization/tdcs'; % Name for the output files
        S.max_total_current = 2e-3; % Maximum total current injected (in mA), optional
        S.max_individual_current = 1e-3; % Maximum current injected per electrode (in mA)
        S.max_active_electrodes = []; % Maximum number of active electrodes in solution, leave empty to unconstrained.
        S.open_in_gmsh = true; % Wether to open output in Gmsh
        S.target = opt_struct('TDCStarget'); % optimization target
        S.avoid = opt_struct('TDCSavoid'); % OPTIONAL: Regions to more strongly avoid

    case 'TDCStarget'
        S.positions = []; % Target positions, in subject space
        S.intensity = 0.2; % Electric field to be reached in target (V/m)
        S.directions = 'normal'; % Electric field direction in target. Cam be either a vector or the string 'normal'
        S.radius = 2; % Radius of target region, in mm
        S.max_angle = []; % Maximum angle between electric field and target direction
        S.indexes = []; % OPTIONAL: Nodes/element indexes of target, overwires "positions"
        S.tissues = []; % OPTINAL: Tissues where the target is defined
        
    case 'TDCSavoid'
        S.positions = []; % Positions to avoid in subject space
        S.weight = 1e3; % Weight for the region. The larger the more it will be avoided
        S.radius = 2; % Radius of avoid region, in mm
        S.indexes = []; % OPTIONAL: Nodes/element indexes of the region, overwires "positions"
        S.tissues = []; % OPTINAL: Tissues where the region is defined

    case 'TDCSDistributedOptimize'
        S.leadfield_hdf = ''; % Name of HDF5 file with leadfield
        S.name = 'optimization/tdcs'; % Name for the output files
        S.max_total_current = 2e-3; % Maximum total current injected (in mA), optional
        S.max_individual_current = 1e-3; % Maximum current injected per electrode (in mA)
        S.max_active_electrodes = []; % Maximum number of active electrodes in solution, leave empty to unconstrained.
        S.target_image = ''; % Image to be "reproduced" via the optimization
        S.mni_space = true; % Wether the image is in MNI space. Set to "false" if in subject space
        S.subpath = ''; % Path to the m2m_{subID} directory. Mandatory if mni_space=true;
        S.intensity = 0.2;  % Target field intensity
        S.min_image_value = 0; %minimum image value to be considered. Corresponds to T_min in Ruffini et al. 2014
        S.open_in_gmsh = true; % Wether to open output in Gmsh

    case 'TmsFlexOptimization'
        S.fnamecoil = ''; % usually to chose from inside resources/coil_models/...
        S.fnamehead = ''; % same as ${subID}.msh created by charm
        S.subpath = ''; % path to the 'm2m_{subID}' folder created by charm (OPTIONAL, filled from fnamehead)
        S.path_optimization = ''; % path to save the results
        S.eeg_cap = ''; % file name of the CSV file with electrode positions; (OPTIONAL, filled from fnamehead)
        S.method = 'distance'; % 'distance' or 'emag'
        S.run_global_optimization = true; % whether to run the DIRECT global optimization
        S.run_local_optimization = true; % whether to run the L-BFGS-B local optimization
        S.run_simulation = true; % whether to run a normal FEM simulation to get the e-field for the final parameters
        S.global_translation_ranges = []; % translation ranges in [mm]; format [[xmin xmax]; [ymin ymax]; [zmin zmax]];
        S.global_rotation_ranges = []; % rotation ranges in [deg]; format [[xmin xmax]; [ymin ymax]; [zmin zmax]];
        S.dither_skip = []; % dithering steps for the grid of the intersection tests (OPTIONAL, standard: 6)
        S.distance = 4; % Distance from the coil to the head in [mm]
        S.fem_evaluation_cutoff = []; % If the penalty from the intersection and self intersection is greater 
                                      % than this cutoff value, the fem will not be evaluated to save time. 
                                      % (OPTIONAL, standard: 1000)
        S.pos = []; % sim_struct.POSITION to define the starting position (optional when a roi is provided)
        S.roi = []; % opt_struct.ROI to define a brain target (optional for method = 'distance')

        % options for DIRECT solver (OPTIONAL)
        % see online documenation of scipy.optimize.direct for a list of available options
        S.direct_args.maxiter = []; % Maximum number of iterations; added just as example here (OPTIONAL, standard 1000)

        % options for L-BFGS-B solver (OPTIONAL)
        % see online documenation of scipy.optimize.minimize(method=’L-BFGS-B’)) for a list of available options
        S.l_bfgs_b_args.options.maxls = []; % maximal line-search steps; added just as example here (OPTIONAL, standard: 100)

    case 'RegionOfInterest'
        S.method = ''; %  method to create the ROI {"manual", "custom", "surface", "volume", "volume_from_surface", "mesh+mask"}
        S.subpath = ''; % path to the 'm2m_{subID}' folder created by charm
        S.mesh = ''; % Path to a mesh
        S.fname_visu = ''; % mesh filename for ROI visualization 

        S.mask_space = ''; %  space the mask is defined in, method = "surface" : {"subject", "subject_lh", "fs_avg_lh", "subject_rh", "fs_avg_rh", "mni"} | method = "volume" : {"subject", "mni"}
        S.mask_path = ''; % path to the mask, method = "surface" : (label, annot, curv, nifti) | method = "volume" : (nifti) (example: "path/to/file")
        S.mask_value = ''; % values that are considered as the mask in the mask files, default 1 (example: 1 | [1, 2])"""
        S.mask_operator = ''; % operator to combine the mask with the ROI {"union", "intersection", "difference"}, default "union"

        S.roi_sphere_center = []; % Nx3 array; sphere center coordinates for spherical masks in mm (example: [0,0,0] | [[1,2,3]; [4,5,6]])
        S.roi_sphere_radius = ''; % Nx1 array; radius of the spherical masks in mm (example: 5 | [3, 45])
        S.roi_sphere_center_space = ''; % space the center coordinates of the spheres are defined in {"subject", "mni"}
        S.roi_sphere_operator = ''; % operator to combine the mask with the ROI {"union", "intersection", "difference"}, default "union"

        S.nodes = []; % Nx3 array; Only for method = "custom" -> a custom list of node coordinates (example: [[1.1, 2.2, 3.3], [4.4, 5.5, 6.6]])"""
        
        S.surface_type = ''; % only for method = "surface" -> Weather the surface is the subject specific central gm surface or a custom surface {"central", "custom"}
        S.surface_path = ''; % only for method = "surface" -> Only for surface_type = "custom" -> The path to a custom surface (msh, freesurfer, gifti)
        S.tissues = []; % Nx1 array of int; only for method = "volume" -> a number of volume tissue tags, default 2 (example: ElementTags.GM | [ElementTags.WM, ElementTags.GM])
        S.surface_inclusion_radius = []; % 1x1 array; only for method = "volume_from_surface" -> The radius from the surface nodes at which the volume elements should be included in mm (example: 5)
        S.node_mask = []; % Nx1 array of bool; only for method = "mesh+mask" -> a boolean node mask (exclusive with elm_mask) (example: [True, ..., False])
        S.elm_mask = []; % Nx1 array of bool; only for method = "mesh+mask" -> a boolean node mask (exclusive with node_mask) (example: [True, ..., False])
        
    case 'TesFlexOptimization'
        S.subpath = ''; % path to the 'm2m_{subID}' folder created by charm
        S.output_folder = ''; % path to save the results
        S.fn_eeg_cap = ''; % EEG cap used for substituting EEG position names with 3d coordinates
        S.goal = ''; % goal function: 'intensity' or 'focality'
        S.threshold = []; % thresholds for focality optimization (vector with length 2)
        S.e_postproc = ''; % postprocessing of e-fields ('magn': magnitude, % 'normal': normal component, 'tangential': tangential component)
        S.electrode = []; % opt_struct.CircularArray or opt_struct.ElectrodeArrayPair
        S.roi = {}; % either a single opt_struct.ROI to define a brain target, or cell array of two opt_struct.ROI {roi, non_roi}
        
    case 'CircularArray' % Nx1 center surround montage                     
        S.radius_inner = []; % radius of inner electrode
        S.radius_outer = []; % radius of outer electrodes
        S.distance_bounds = []; % distance bounds (min, max) between inner and outer electrodes
        S.n_outer = []; % number of outer electrodes
        S.dirichlet_correction = ''; % set to true when all outer electrodes are connected to the same channel (slower)
        S.current = []; % currents for each electrode

   case 'ElectrodeArrayPair' % Pair of TES electrode arrays
        S.length_x = []; % x-dimension of electrodes
        S.length_y = []; % y-dimension of electrodes
        S.dirichlet_correction_detailed = false;  % account for inhomogenous current distribution at electrode-skin interface (slow)
        S.current = [];  % currents for each electrode      
        
end
