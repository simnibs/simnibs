%
% Example script that shows how to extract results from
% cortical regions-of-interest using node_data indices.
% 
% It asks for simulation results saved in the fsavg_overlays
% subfolder of a simulation, and reports some key values
%       (1) for the whole cortex
%       (2) for the (randomly selected) "S_precentral-sup-part" area 
%            (from the Destrieux atlas or "a2009s" in FreeSurfer notation) 
%            of the left hemisphere
%       (3) in a sphere of 10 mm radius around the peak position
%
% A. Thielscher, 06-Nov-2017; updated 08-Oct-2018
% G. Saturnino, updated 18-SEP-2019


% get name of results mapped on fsaverage surface
[fname,pname] = uigetfile('*_fsavg.msh', ...
                          'Select a results file in a fsavg_overlays folder');
if isequal(fname,0) || isequal(pname,0); return; end

% Load simulations
m=mesh_load_gmsh4(fullfile(pname,fname));

% Show all the fields
for i =1:length(m.node_data)
    mesh_show_surface(m, 'field_idx',i)
end

% load the fsaverage surface with the Destrieux atlas:
% 	labels is a mesh structure containing the surface and label data
%   s is the list of label names
[labels, snames]=mesh_load_fssurf('fsaverage','label','a2009s');


% show the atlas (loaded as first --> node data index 1)
mesh_show_surface(labels,'field_idx',1)


% -------------------------------------------------------
% EXAMPLE 1 
% display some key results for whole cortex
% Note: Focality results are only an estimate, as the fsaverage surface
% rather the individual surface was used to determine the stimulated areas
disp(' ')
disp('whole cortex:')
% You can modify the 'field_idx' to see different fields
summary=mesh_get_fieldpeaks_and_focality(m,'field_idx',1);


% -------------------------------------------------------
% EXAMPLE 2
% extract "S_precentral-sup-part" of the left hemisphere

% find the atlas index of the area
area_idx=find(strcmpi(snames,'lh.S_precentral-sup-part'));
% find all nodes having this atlas index
node_idx=labels.node_data{1}.data==area_idx;
% extract those nodes and the related triangles and data
m_ROI=mesh_extract_regions(m, 'node_idx', node_idx);

% show the extracted ROI:
% show the whole cortex semi-transparent
mesh_show_surface(labels,'showSurface',true,'facealpha',0.3); 
% add the extracted area to the plot
mesh_show_surface(m_ROI,'showSurface',true,'surfaceColor',[1 0 0],'haxis',gca); 
title('lh.S_precentral-sup-part');

% get some key results for the "S_precentral-sup-part" area
disp(' ')
disp('lh.S_precentral-sup-part:')
mesh_get_fieldpeaks_and_focality(m_ROI,'field_idx',1);


% -------------------------------------------------------
% EXAMPLE 3
% extract a spherical ROI with 10 mm radius around the peak position
%
% the "summary" structure contains the peak values
% together with their positions:
%   "summary.percentiles" lists the tested percentile cutoffs - 
%   the 99.9 percentile is the 3rd entry
%
%   "summary.perc_values" lists the corresponding values, as they
%   are also displayed by mesh_get_fieldpeaks_and_focality
%
%   "summary.XYZ_perc" lists the corresponding center positions -
%   the center position for the 99.9 percentile is the 3rd row:
peak_pos=summary.XYZ_perc(3,:);

% distance to peak position
dist=sqrt(sum(bsxfun(@minus,m.nodes,peak_pos).^2,2));

% extract nodes closer than 10 mm, and the related triangles and data
node_idx=dist<10;
m_ROI=mesh_extract_regions(m, 'node_idx', node_idx);

% show the extracted ROI:
% show the whole cortex semi-transparent
mesh_show_surface(m,'showSurface',true,'facealpha',0.3); 
% add the extracted area to the plot
mesh_show_surface(m_ROI,'showSurface',true,'surfaceColor',[1 0 0],'haxis',gca); 
title('10 mm spherical ROI around peak position');

% get some key results for the spherical ROI
disp(' ')
disp('10 mm spherical ROI around peak position:')
mesh_get_fieldpeaks_and_focality(m_ROI,'field_idx',1);



