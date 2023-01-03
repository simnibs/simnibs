function varargout = mesh_show_surface(m, varargin)
% plots a surface or results on a surface
%
% USAGE:
%   hfig = mesh_show_surface(m [,'OptionName',OptionValue,...])
%
%   m: input mesh
%
%   Options are set using the option name followed by the value.
%   Supported options:
%   field_idx: identifier of a results field;
%              can be either the name or the index of the
%              element_data or node_data field (standard: 'magnE', or the
%              first in case only one data field is loaded)
%   region_idx: uses the triangle surface(s) with the given region numbers
%               (standard: 1002 for gray matter for .msh-files or all regions
%                for Freesurfer- and CAT12-surface-files)
%   datatype: 'node' for node_data or 'tri' for triangle element_data
%                (standard: 'tri', except element_data is empty)
%   showElec: show electrodes as transparent surfaces (standard: true)
%   showSurface: show only surface instead of data (standard: false)
%   surfaceColor: surface color as RGB triplet (standard: [0.8 0.8 0.8])
%   scaleLimits: upper and lower limits of the colorbar
%               (standard: [0 absmax] for positive-only data;
%                          [-absmax absmax] for other data
%                absmax is max(abs(0.1 percentile),99.9 percentile)))
%   colormap: colormap for data (standard: jet for positive-only data;
%                                blue-gray-red for other data
%   facealpha: transparency setting (from 0 -fully tranparent- to 1 -fully
%   opaque-; standard: 1)
%   haxis: handle to an existing axis, which will be used for plotting
%
%   Output:
%   hfig: handle of the created figure (optional)
%
% Examples:
%  mesh_show_surface(m); % plots magnE on the gray matter surface
%  mesh_show_surface(m,'field_idx','magnJ','region_numbers',1001); % plots magnJ on white matter
%  mesh_show_surface(m,'showElec',false); % doesn't show the electrodes
%         
% A. Thielscher 09-Sep-2018

% standard settings and behavior
s.field_idx = 'magnE';
s.region_idx = 1002;
s.datatype = 'tri';
s.showElec = true;
s.showSurface = false;
s.surfaceColor = [0.8 0.8 0.8];
s.scaleLimits = [];
s.colormap = [];
s.haxis =[];
s.facealpha=1;

if isempty(m.element_data); s.datatype = 'node'; end
if isempty(m.tetrahedra); s.region_idx=unique(m.triangle_regions)'; end % work-around to guess that input was not a .msh file
if strcmpi(s.datatype,'node')
    if isempty(m.node_data); s.showSurface = true; end
    if length(m.node_data)==1; s.field_idx=1; end
else
    if isempty(m.element_data); s.showSurface = true; end
    if length(m.element_data)==1; s.field_idx=1; end
end

% parse input
if nargin<1; error('mesh is needed as input'); end
s=parse_input(s,varargin{:});

% get index of data field in case it was given as field name
if ~s.showSurface; s.field_idx = get_field_idx(m,s.field_idx,s.datatype); end

% extract region
disp(['Using region number(s) ' num2str(s.region_idx)]);
morg=m; % copy needed for showElec
m=mesh_extract_regions(m,'elemtype','tri','region_idx',s.region_idx);

% get data, name and lower and upper limits
if ~s.showSurface
    [data, name, s.scaleLimits] = get_data_and_scaleLimits(m,s.field_idx,s.datatype,s.scaleLimits);
    % scaling is from 0 to 99.9 percentile (magnE, magnJ)
    % or from .1 to 99.9 percentile (normal components)
    % scaling will be only updated when empty
    
    if isempty(s.colormap)
        if s.scaleLimits(1)>=0; s.colormap='jet';
        else
           cmap=jet(64);
           cmap=[cmap(24:-1:1,:);zeros(5,3);cmap(end:-1:41,:)];
           weight_cmap = ones(53,1);
           weight_cmap(27-8:27)=(8:-1:0)/9;
           weight_cmap(27:27+8)=(0:8)/9;
           for i=1:3
               cmap(:,i)=weight_cmap.*cmap(:,i) + (1-weight_cmap).*0.7;
           end
           s.colormap=cmap;
        end
    end
end


% plot
% ----
if isempty(s.haxis); hfig=figure;
else; subplot(s.haxis); end

if s.showSurface
    hp=patch('Faces',m.triangles,'Vertices',m.nodes,'FaceVertexCData',s.surfaceColor,...
          'FaceColor','flat','EdgeColor','none','FaceAlpha',s.facealpha);
    axis equal
    axis off
else
    if min(size(data))>1; error('Scalar data required'); end
    
    if strcmpi(s.datatype,'node')
        hp=patch('Faces',m.triangles,'Vertices',m.nodes,'FaceVertexCData',data,...
                 'FaceColor','interp','EdgeColor','none','CDataMapping','scaled','FaceAlpha',s.facealpha);
    else
        hp=patch('Faces',m.triangles,'Vertices',m.nodes,'FaceVertexCData',data,...
                 'FaceColor','flat','EdgeColor','none','CDataMapping','scaled','FaceAlpha',s.facealpha);
    end
    axis equal
    axis off
    colormap(s.colormap);
    h=colorbar;
    h=get(h,'label');
    set(gca,'CLim',s.scaleLimits);
    
    title(name,'Interpreter','none');        
    if strcmpi(name,'magnE')
        set(h,'String','electric field strength in [V/m]');
    elseif strcmpi(name,'magnJ')
        set(h,'String','current density in [A/m²]');
    elseif strcmpi(name,'E.normal')
        set(h,'String','normal component of electric field in [V/m]');
    elseif strcmpi(name,'J.normal')
        set(h,'String','normal component of current density in [A/m²]');
    end    
end

material(hp,'dull');
lighting gouraud
hlight=camlight('headlight');
set(gca,'UserData',hlight);
hrot = rotate3d;
set(hrot,'ActionPostCallback',@(~,~)camlight(get(gca,'UserData'),'headlight'));

if s.showElec
    hold on
    idxElec=morg.triangle_regions>1500&morg.triangle_regions<2000;
    patch('Faces',morg.triangles(idxElec,:),'Vertices',morg.nodes,'FaceColor','k','FaceAlpha',0.2,'EdgeColor','none');
end

if nargout; varargout{1}=hfig; end
