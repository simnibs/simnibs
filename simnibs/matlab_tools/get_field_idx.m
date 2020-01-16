function field_idx = get_field_idx(m,field_idx,datatype)
% returns index of data field
%
% A. Thielscher, 23-Sep-2018

if strcmpi(datatype,'node'); data_fields=m.node_data;
else; data_fields=m.element_data; end

if ~isnumeric(field_idx)
    for i=1:length(data_fields)
        if strcmpi(data_fields{i}.name,field_idx)
            field_idx=i;
            break;
        end
    end
end

if ~isnumeric(field_idx)
     error(['Data field ' field_idx ' not found in mesh data']);
else
     disp(['Extracting data for ' data_fields{field_idx}.name]);
end

