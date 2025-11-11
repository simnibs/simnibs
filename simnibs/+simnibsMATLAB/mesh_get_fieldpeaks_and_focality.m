function varargout = mesh_get_fieldpeaks_and_focality(m, varargin)
% extracts some key numbers on peak field strength and the focality
% of stimulation
%
% USAGE:
%   s = mesh_get_fieldpeaks_and_focality(m [,'OptionName',OptionValue,...])
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
%   percentiles: used to determine the peak values, given as row vector;
%                for negative data, 100 minus the given values is used
%               (standard: [95 99 99.9] percentiles) 
%   focality_cutoffs: cutoffs in % of 99.9 percentile; the volume or area
%                     in which this cutoff is exceeded is reported as index
%                     of focality (standard:[50 75 100]);
%                     for negative data, cutoffs apply to the 0.1 percentile
%   printsummary: print summary of results (standard: true)
%
% OUTPUT:
%  s is a structure with the following entries:
%   field_idx: index of the used element_data or node_data field
%   field_name: name of the data field
%   region_idx: indices of the used regions
%   datatype: 'node', 'tri' or 'tet'
%   sizeunits: units used for the focality results; either 'square mm' 
%              for surface data or 'cubic mm' for volume data
%   valueunits: units of the extracted data; either V/m, A/m� or empty
%               (the latter for unknown data)
%   max: maximal value; this value can be caused by an outlier; please
%        visually check that it is embedded in a region of high values
%        in case you want to use it
%   percentiles: list of used percentiles, and
%   perc_values: corresponding values
%   focality_cutoffs: cutoffs in % of 99.9 percentile for determining the 
%                     focality, and
%   focality_values: corresponding values
%   XYZ_max: position of the maximum
%   XYZ_perc: mean positions of the elements or nodes contributing to the 
%             listed perc_values
%   XYZstd_perc: standard deviation of the elements or nodes contributing
%                to the listed perc_values
%
%  In case the data contains negative values, also the following fields
%  will appear:
%   min: minimal value; see remarks for max regarding its robustness
%   perc_neg_values: values for the lower percentiles of the distribution
%                    (standard: for the [5 1 0.1] percentiles).
%   focality_neg_values: focality data for the negative values
%   XYZ_min: position of the minimum
%   XYZ_perc_neg: mean positions of the elements or nodes contributing
%                 to the listed perc_neg_values
%   XYZstd_perc_neg: standard deviation of the elements or nodes
%                    contributing to the listed perc_neg_values    
%
% Examples:
%  s=mesh_get_fieldpeaks_and_focality(m); % plot summary for GM, return results structure
%  mesh_get_fieldpeaks_and_focality(m,'region_idx',1); plot summary for WM (.msh files)
%  s=mesh_get_fieldpeaks_and_focality(m'printsummary',false); % do not plot summary
%         
% A. Thielscher 07-Sep-2018

% standard settings and behavior
s.field_idx = 'magnE';
s.region_idx = 2;
s.datatype = 'tet';
s.percentiles=[95 99 99.9];
s.focality_cutoffs=[50 75];
s.printsummary=true;

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
s.sizeunits='cubic mm';
if strcmpi(s.datatype,'tet')
    m=mesh_extract_regions(m,'elemtype','tet','region_idx',s.region_idx);
elseif strcmpi(s.datatype,'tri')
    m=mesh_extract_regions(m,'elemtype','tri','region_idx',s.region_idx);
    s.sizeunits='square mm';
else
    m=mesh_extract_regions(m,'elemtype','both','region_idx',s.region_idx);
    if isempty(m.tetrahedra); s.sizeunits='square mm'; end
end

% get data, field name, scaleLimits, element sizes and positions
[data, name, scaleLimits, elemsizes, elempos] = get_data_and_scaleLimits(m,s.field_idx,s.datatype,[]);

if min(size(data))>1; error('Scalar data required'); end
idx=~isnan(data); % NaNs occured for .angle data; make sure to get rid of them
data=data(idx);
elemsizes=elemsizes(idx);
elempos=elempos(idx,:);

% sort data, get cdf
[data,idx] = sort(data);
elemsizes=elemsizes(idx);
elemsizes=cumsum(elemsizes);
elemsizesNormed=elemsizes/elemsizes(end);
elempos=elempos(idx,:);


% extract peak and percentile values
s.field_name=name;
s.valueunits='';
if strcmpi(name,'magnE')||strcmpi(name,'E.normal')
    s.valueunits='[V/m]';
elseif strcmpi(name,'magnJ')||strcmpi(name,'J.normal')
    s.valueunits='[A/m�]';
end
s.max=max(data);
s.XYZ_max=elempos(data==s.max,:);

for i=1:length(s.percentiles)
    idx=find(elemsizesNormed>s.percentiles(i)/100,1,'first');
    if isempty(idx); idx=length(data); end
    s.perc_values(i)=data(idx);
    
    % mean and SD of positions (equations weighted for element size)
    meanVal=sum(bsxfun(@times,elempos(idx:end,:),elemsizes(idx:end)),1)./...
            repmat(sum(elemsizes(idx:end)),1,3);
    
    N_nonzero=sum(elemsizes(idx:end,:)>0);
    scaleFac=(N_nonzero-1)/N_nonzero*sum(elemsizes(idx:end,:));
    for j=1:3
        s.XYZstd_perc(i,j)=sqrt( sum(elemsizes(idx:end).*( elempos(idx:end,j)-meanVal(j) ).^2,1)/scaleFac );
    end
    s.XYZ_perc(i,:)=meanVal;    
end

if scaleLimits(1)<0
    s.min=min(data);
    s.XYZ_min=elempos(data==s.min,:);
    
    perc_neg=100-s.percentiles;
    for i=1:length(perc_neg)
        idx=find(elemsizesNormed<perc_neg(i)/100,1,'last');
        if isempty(idx); idx=1; end
        s.perc_neg_values(i)=data(idx);
                
        % mean and SD of positions (equations weighted for element size)
        meanVal=sum(bsxfun(@times,elempos(1:idx,:),elemsizes(1:idx)),1)./...
                repmat(sum(elemsizes(1:idx)),1,3);
    
        N_nonzero=sum(elemsizes(1:idx,:)>0);
        scaleFac=(N_nonzero-1)/N_nonzero*sum(elemsizes(1:idx,:));
        for j=1:3
            s.XYZstd_perc_neg(i,j)=sqrt( sum(elemsizes(1:idx).*( elempos(1:idx,j)-meanVal(j) ).^2,1)/scaleFac );
        end
        s.XYZ_perc_neg(i,:)=meanVal;
    end
end


% extract focality
idx=find(elemsizesNormed>99.9/100,1,'first');
if isempty(idx); idx=length(data); end 
peakvalue=data(idx);
for i=1:length(s.focality_cutoffs)
    idx=find(data>=s.focality_cutoffs(i)/100*peakvalue,1,'first');
    if idx == 1
        s.focality_values(i) = elemsizes(end);
    else
        s.focality_values(i)=elemsizes(end)-elemsizes(idx-1);
    end
end

if scaleLimits(1)<0
    idx=find(elemsizesNormed<0.1/100,1,'last');
    if isempty(idx); idx=1; end 
    peakvalue=data(idx);
    for i=1:length(s.focality_cutoffs)
        idx=find(data<=s.focality_cutoffs(i)/100*peakvalue,1,'last');
        s.focality_neg_values(i)=elemsizes(idx);
    end    
end


% print summary
if s.printsummary
    disp(' ');
    disp('---------------------------------------------');
    disp('SUMMARY ');
    disp(['field name: ',s.field_name]);
    disp(['region indices: ',num2str(s.region_idx)]);
    disp(' ');
    disp('peak fields ');
    disp(['percentiles:      ',num2str(s.percentiles,3)]);
    disp(['values:           ',num2str(s.perc_values,3),' (in ',s.valueunits,')']);
    if scaleLimits(1)<0
        disp(['values (negative): ',num2str(s.perc_neg_values,3)]);
    end
    disp(' ');
    disp('focality');
    disp(['cutoffs:          ',num2str(s.focality_cutoffs,3), ' (in % of 99.9 percentile)']);
    disp(['values:           ',num2str(s.focality_values,3), ' (in ', s.sizeunits, ')']);
    if scaleLimits(1)<0
        disp(['values (negative): ',num2str(s.focality_neg_values,3)]);
    end
    disp('---------------------------------------------');
end

% prepare output
if nargout
    s=rmfield(s,'printsummary');
    if scaleLimits(1)<0
        field_str={'field_idx','field_name','region_idx','datatype','sizeunits','valueunits',...
                   'max','percentiles','perc_values','focality_cutoffs','focality_values',...
                   'min','perc_neg_values','focality_neg_values',...
                   'XYZ_max','XYZ_perc','XYZstd_perc',...
                   'XYZ_min','XYZ_perc_neg','XYZstd_perc_neg'};
    else
        field_str={'field_idx','field_name','region_idx','datatype','sizeunits','valueunits',...
                   'max','percentiles','perc_values','focality_cutoffs','focality_values',...
                   'XYZ_max','XYZ_perc','XYZstd_perc'};
    end
    s=orderfields(s,field_str);
    varargout{1}=s;
end

% % for testing
% figure
% h1 = subplot(1,2,1);
% mesh_show_surface(m,'haxis',h1)
% h2= subplot(1,2,2);
% mesh_show_surface(m,'haxis',h2,'showSurface',true,'facealpha',0.1)
% hold on
% [X,Y,Z]=ellipsoid(s.XYZ_max(1),s.XYZ_max(2),s.XYZ_max(3),2,2,2);
% surf(X,Y,Z,'LineStyle','none','FaceColor',[1 0 0])
% for i=1:size(s.XYZ_perc,1)
%     [X,Y,Z]=ellipsoid(s.XYZ_perc(i,1),s.XYZ_perc(i,2),s.XYZ_perc(i,3),...
%                       s.XYZstd_perc(i,1),s.XYZstd_perc(i,2),s.XYZstd_perc(i,3));
%     surf(X,Y,Z,'LineStyle','none','FaceColor',[1 0 0],'facealpha',0.2)
% end
% if scaleLimits(1)<0
%     [X,Y,Z]=ellipsoid(s.XYZ_min(1),s.XYZ_min(2),s.XYZ_min(3),2,2,2);
%     surf(X,Y,Z,'LineStyle','none','FaceColor',[0 0 1])
%     for i=1:size(s.XYZ_perc_neg,1)
%         [X,Y,Z]=ellipsoid(s.XYZ_perc_neg(i,1),s.XYZ_perc_neg(i,2),s.XYZ_perc_neg(i,3),...
%                           s.XYZstd_perc_neg(i,1),s.XYZstd_perc_neg(i,2),s.XYZstd_perc_neg(i,3));
%         surf(X,Y,Z,'LineStyle','none','FaceColor',[0 0 1],'facealpha',0.2)
%     end
% end
