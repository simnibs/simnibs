function varargout = mesh_get_simulation_result
% 
% The function plots the electric field distribution on GM (field strength,
% or normal component) and determines some basic key values on peak
% stimulation strength and focality.
% 
% It reads a SimNIBS mesh with simulation results or the results in the
% subject_overlays or fsavg_overlays subfolders.
% 
% * It plots a histogram of the field values and the field on the
%   GM surface.
%
% * To characterize the stimulation intensity, it determines the 
%   values for the 95, 99 and 99.9 percentiles.
%
% * In case the results contain negative values (e.g. for the normal
%   component), it also determines the values for the 5, 1 and 0.1
%   percentiles.
% 
% * To characterize focality, it determines the volume (for .msh files) or
%   surface area (for the results in subject_overlays or fsavg_overlays) in 
%   which the field exceeds certain fractions of the 99.9 percentile 
%   (fractions of the 0.1 percentile for negative values).
%  (Please note that focality values for results in fsavg_overlays
%   are based on the area of the fsaverage template. It is thus
%   preferable to report focality results for .msh files or for results
%   in the subject_overlays folder)
%
% In order to adapt it to your needs, please see the help of the underlying
% functions mesh_get_histogram, mesh_show_surface and 
% mesh_get_fieldpeaks_and_focality
% 
% A. Thielscher, 11-Apr-2018; updated 13-Sep-2018


% get name of mesh or of results mapped on surface
[fname,pname] = uigetfile('*.msh;lh.*.magn;lh.*.normal;lh.*.tangent;lh.*.angle', ...
                          'Select .msh file or results in subject_overlays or fsavg_overlays');
if isequal(fname,0) || isequal(pname,0); return; end


% load
[~,~,ext] = fileparts(fname);
if strcmpi(ext,'.msh')
    m=mesh_load_gmsh4(fullfile(pname,fname));
else
    m=mesh_load_fsresults(fullfile(pname,fname));
end
    

% plot histogram of field in GM (based on volume for .msh; area for surface data) 
mesh_get_histogram(m)

% show field on GM surface
mesh_show_surface(m)

% extract a few key values
mesh_get_fieldpeaks_and_focality(m);

if nargout; varargout{1}=s; end
    


