%
% Example script that demonstrates how to access data in .msh-files
% 
% It reads a SimNIBS mesh with simulation results and reports some key values
%       (1) for the whole cortex
%       (2) in a sphere of 10 mm radius around the peak position
%
% It plots the field on the GM surface and the extracted region of interest
%
% A. Thielscher, 11-Apr-2018; updated 13-Sep-2018


% load mesh
[fname,pname] = uigetfile('*.msh','Select a mesh with simulation results');
if isequal(fname,0) || isequal(pname,0); return; end
m=mesh_load_gmsh4([pname fname]);


% -------------------------------------------------------
% EXAMPLE 1 
% display some key results for whole cortex
disp(' ')
disp('whole cortex:')
summary=mesh_get_fieldpeaks_and_focality(m,'field_idx',2);

% show field on the GM surface
mesh_show_surface(m)


% -------------------------------------------------------
% EXAMPLE 2
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

% extract nodes closer than 10 mm, and belonging to tetrahedra with region
% number 2 (gray matter has region number 2)
node_idx=dist<10;
m_ROI=mesh_extract_regions(m, 'elemtype','tet','region_idx',2,'node_idx', node_idx);

% get some key results for the spherical ROI
disp(' ')
disp('10 mm spherical ROI around peak position:')
mesh_get_fieldpeaks_and_focality(m_ROI,'field_idx',2);

% show the extracted ROI:
% show the whole cortex semi-transparent
mesh_show_surface(m,'showSurface',true,'facealpha',0.4);
% create new mesh containing the enclosing surface of the ROI tetrahedra
TR = triangulation(m_ROI.tetrahedra,m_ROI.nodes(:,1),m_ROI.nodes(:,2),m_ROI.nodes(:,3));
m_vis=mesh_empty;
m_vis.nodes=m_ROI.nodes;
m_vis.triangles=freeBoundary(TR);
m_vis.triangle_regions=ones(size(m_vis.triangles,1),1);
% add the he enclosing surface to the plot
mesh_show_surface(m_vis,'showSurface',true,'surfaceColor',[1 0 0],'haxis',gca); 
title('10 mm spherical ROI around peak position');