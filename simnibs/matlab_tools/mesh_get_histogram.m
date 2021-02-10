function varargout = mesh_get_histogram(m, varargin)
% creates a histogram of the field distribution
%
% USAGE:
%   [BinCenters BinData] = mesh_get_histogram(m [,'OptionName',OptionValue,...])
%
%   m: input mesh
%
%   Options are set using the option name followed by the value.
%   Supported options:
%   field_idx: identifier of a results field;
%              can be either the name or the index of the
%              element_data or node_data field (standard: 'magnE', or the
%              first in case only one data field is loaded)
%   region_idx: uses the elements with the given region numbers
%               (standard: 2 for gray matter tetrahedra or all triangles
%                for Freesurfer- and CAT12-surface-files)
%   datatype: 'node' for node_data, 'tri' for triangle element_data,
%              'tet' for tet element_data
%              (standard: 'tet', except element_data is empty)
%   relscale: scale histogram relative to the volume size
%             (optional; standard: false; output will then be in mm³)
%   scaleLimits: upper and lower limits of the histogram
%               (standard: [0 absmax] for positive-only data;
%                          [-absmax absmax] for other data
%                absmax is max(abs(0.1 percentile),99.9 percentile)))
%   nbins: number of bins (optional; standard: 100 bins)
%   haxis: handle to an existing axis, which will be used for plotting
%
%   Output:
%   BinCenters: bin positions on the x axis
%   BinData: the histogram data
%   (BinData and BinCenters are optional; the histogram will be plotted 
%    in a figure when they are not defined)
%
% Examples:
%  mesh_get_histogram(m); % plot histogram for magnE for gray matter
%  [BinCenters BinData] =  mesh_get_histogram(m,'magnJ',1); % return histogram data 
%                                                           % for magnJ in white matter
%         
% A. Thielscher 07-Sep-2018

% standard settings and behavior
s.field_idx = 'magnE';
s.region_idx = 2;
s.relscale = false;
s.scaleLimits = [];
s.datatype = 'tet';
s.nbins = 100;
s.haxis = [];

if isempty(m.element_data); s.datatype = 'node'; end
if isempty(m.tetrahedra); s.region_idx=unique(m.triangle_regions)'; end % work-around to guess that input was not a .msh file

if strcmpi(s.datatype,'node')
    if length(m.node_data)==1; s.field_idx=1; end
else
    if length(m.element_data)==1; s.field_idx=1; end
end

% parse input
if nargin<1; error('mesh is needed as input'); end
s=parse_input(s,varargin{:});

% get index of data field in case it was given as field name
s.field_idx = get_field_idx(m,s.field_idx,s.datatype);

% extract regions
disp(['Using region number ' num2str(s.region_idx)]);
if strcmpi(s.datatype,'tet')
    m=mesh_extract_regions(m,'elemtype','tet','region_idx',s.region_idx);
elseif strcmpi(s.datatype,'tri')
    m=mesh_extract_regions(m,'elemtype','tri','region_idx',s.region_idx);
else
    m=mesh_extract_regions(m,'elemtype','both','region_idx',s.region_idx);
end

% get data, field name, scaleLimits (when they are empty) and element sizes
[data, name, s.scaleLimits, elemsizes] = get_data_and_scaleLimits(m,s.field_idx,s.datatype,s.scaleLimits);
% scaling is from 0 to 99.9 percentile (magnE, magnJ)
% or from .1 to 99.9 percentile (normal components)
% scaling will be only updated when empty


% create histogram
% --------------------------
if min(size(data))>1; error('Scalar data required'); end

BinWidth=(s.scaleLimits(2)-s.scaleLimits(1))/s.nbins;
LowerEdges=s.scaleLimits(1) + [0:s.nbins-1]*BinWidth;
BinCenters=LowerEdges+BinWidth/2;

BinData=zeros(s.nbins,1);
BinData(1) = sum(elemsizes(data<LowerEdges(1)+BinWidth));
for i=2:length(LowerEdges)-1
    BinData(i) = sum(elemsizes((data>=LowerEdges(i))&(data<LowerEdges(i)+BinWidth)));
end
BinData(end) = sum(elemsizes(data>=LowerEdges(end)));
if s.relscale; BinData=BinData./sum(elemsizes); end

if nargout
    varargout{1}=BinCenters;
    varargout{2}=BinData;
else
    if isempty(s.haxis); figure;
    else; subplot(s.haxis); end
    
    bar(BinCenters,BinData,1,'edgecolor','none');
    title(['pdf of ' name]);
    
    % y axis labeling
    if s.relscale; ylabel('frequency of occurence');
    elseif strcmpi(s.datatype,'tet'); ylabel('volume in [mm³]'); 
    elseif strcmpi(s.datatype,'tri'); ylabel('area in [mm²]'); 
    elseif strcmpi(s.datatype,'node')
        if isempty(m.triangles); ylabel('volume in [mm³]');
        else;  ylabel('area in [mm²]'); end
    end
    
    % x axis labeling
    if strcmpi(name,'magnE')
        xlabel('electric field strength in [V/m]');
    elseif strcmpi(name,'magnJ')
        xlabel('current density in [A/m²]');
    elseif strcmpi(name,'E.normal')
        xlabel('normal component of electric field in [V/m]');
    elseif strcmpi(name,'J.normal')
        xlabel('normal component of current density in [A/m²]');
    end
end
    
    
   