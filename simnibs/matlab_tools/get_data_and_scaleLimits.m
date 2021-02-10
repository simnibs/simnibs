function [data, name, scaleLimits, varargout] = get_data_and_scaleLimits(m,field_idx,datatype,scaleLimits)
% returns the data, the field name and fills in the scaleLimits (if they
% are empty); optionally, also the element sizes are returned
% 
% USAGE:
%  [data, name, scaleLimits] = get_data_and_scaleLimits(m,field_idx,datatype,scaleLimits)
%  [data, name, scaleLimits, elemsizes] = get_data_and_scaleLimits(m,field_idx,datatype,scaleLimits)
%
% The scale limits are determined as follows:
% [0 absmax] for positive-only data;
% [-absmax absmax] for other data
% absmax is max(abs(0.1 percentile),99.9 percentile)))
%
% NOTE: To get scaleLimits or elemement sizes for node_data,
% the mesh must contain either tetrahedra or triangles. Otherwise,
% the mapping of sizes (tet volumes, tri areas) to node is ambiguous.
%
% A. Thielscher, 23-Sep-2018

% get data and name
if strcmpi(datatype,'tet')
    data=m.element_data{field_idx}.tetdata;
    name=m.element_data{field_idx}.name;
elseif strcmpi(datatype,'tri')
    data=m.element_data{field_idx}.tridata;
    name=m.element_data{field_idx}.name;
else
    data=m.node_data{field_idx}.data;
    name=m.node_data{field_idx}.name;
end

% get element sizes if needed
if isempty(scaleLimits)||(nargout>3)
    if strcmpi(datatype,'tet')
        % tet volumes
        elemsizes=mesh_get_tetrahedron_sizes(m);
    elseif strcmpi(datatype,'tri')
        % triangle areas
        elemsizes=mesh_get_triangle_sizes(m);
    else % node data
        if ~isempty(m.triangles)&&~isempty(m.tetrahedra)
            disp('Both tets and triangles found: Ambiguity to assign sizes to nodes')
            error('mesh has to contain either tets or triangles, not both');
        end
        if isempty(m.triangles)
            % use tet volumes to get node sizes
            elemsizes=mesh_get_tetrahedron_sizes(m);
            elemsizes=elemsizes/4;
            nodesizes=zeros(size(m.nodes,1),1);
            for i=1:size(m.tetrahedra,1)
                nodesizes(m.tetrahedra(i,:))=nodesizes(m.tetrahedra(i,:))+elemsizes(i);
            end
            elemsizes=nodesizes;
        else
            % use tri areas to get node sizes
            elemsizes=mesh_get_triangle_sizes(m);          
            elemsizes=elemsizes/3;
            nodesizes=zeros(size(m.nodes,1),1);
            for i=1:size(m.triangles,1)
                nodesizes(m.triangles(i,:))=nodesizes(m.triangles(i,:))+elemsizes(i);
            end
            elemsizes=nodesizes;
        end
    end
end

% get element sizes if needed
if nargout>3; varargout{1}=elemsizes; end

% get element postions if needed
if nargout>4
    if strcmpi(datatype,'tet')
        varargout{2}=mesh_get_tetrahedron_centers(m);
    elseif strcmpi(datatype,'tri')
        varargout{2}=mesh_get_triangle_centers(m);
    else
        varargout{2}=m.nodes;
    end
end

% get scaleLimits if needed
if isempty(scaleLimits)
    if size(data,2)>1
        disp('vector data: Using magnitude to get scales');
        dataHlp=sqrt(sum(data.^2,2));
    else
        dataHlp=data;
    end
    
    idx=~isnan(dataHlp);
    dataHlp=dataHlp(idx);
    elemsizesHlp=elemsizes(idx);
    
    [dataHlp,idx] = sort(dataHlp);
    elemsizesHlp=elemsizesHlp(idx);
    elemsizesHlp=cumsum(elemsizesHlp);
    elemsizesHlp=elemsizesHlp/elemsizesHlp(end);
    idx=find(elemsizesHlp>0.999,1,'first');
    if isempty(idx); idx=length(dataHlp); end
    scaleLimits=[0 dataHlp(idx)]; % scale from 0 to 99.9 percentile

    if dataHlp(1)<0
        % rescale from 0.1 to 99.9 percentile
        idx=find(elemsizesHlp<0.001,1,'last');
        if isempty(idx); idx=1; end
        absmax=max([-dataHlp(idx) scaleLimits(2)]);
        scaleLimits = [-absmax absmax];
    end
end