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


validtypes={'TMSoptimize', 'TDCSoptimize', 'TDCStarget', 'TDCSavoid', 'TDCSDistributedOptimize'};

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
       S.fnamecoil =  '';      % to chose from inside ccd-files
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
       S.method = 'direct'; % Solution method, either 'direct' or 'ADM'. The former is only valid with .ccd coil files

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

end

