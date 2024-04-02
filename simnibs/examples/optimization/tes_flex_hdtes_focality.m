% Example to run TESoptimize with an HDTES montage to optimize the field focality in the ROI vs non-ROI in MATLAB
%
% Written by: Konstantin Weise (2023)

% Initialize structure
opt = opt_struct('TESoptimize');

% path of m2m folder containing the headmodel
opt.subpath = "m2m_ernie";

% output folder
opt.output_folder = "tes_optimze_hdtes_focality";

% type of goal function
opt.goal = "focality";

% define threshold(s)
opt.threshold = [0.1, 0.2];

% postprocessing of e-fields ("magn": magnitude, "normal": normal component, "tangential": tangential component)
opt.e_postproc = "magn";

% define electrode
electrode = opt.add_electrode();
electrode.type = "CircularArray";                                     % HDTES center surround montage
electrode.radius_inner = 10;                                          % radius of inner electrode
electrode.radius_outer = 10;                                          % radius of outer electrodes
electrode.distance_bounds = [25, 100];                                % distance bounds between inner and outer electrodes
electrode.n_outer = 4;                                                % number of outer electrodes
electrode.dirichlet_correction = false;                               % electrode wise dirichlet correction
electrode.dirichlet_correction_detailed = false;                      % node wise dirichlet correction
electrode.current = [0.002, -0.002/4, -0.002/4, -0.002/4, -0.002/4];  % initial currents

% define ROI
roi = opt.add_roi();
roi.type = "GMmidlayer";

% center of spherical ROI in subject space (in mm)
roi.roi_sphere_center_subject = [-41, -13,  66];

% radius of spherical ROI (in mm)
roi.roi_sphere_radius = 20;

% define non-ROI
roi = opt.add_roi();
roi.type = "custom";
roi.domains = [1, 2];

% Run optimization
run_simnibs(opt);
