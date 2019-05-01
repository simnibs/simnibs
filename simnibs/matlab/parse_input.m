function s=parse_input(s,varargin)
% parses the input arguments into the structure s
%
% A. Thielscher, 10-Sep-2018

names = fieldnames(s);

if ~mod(nargin,2)
    error('options have to be given as pairs of option name and option value');
end

for i=1:2:nargin-1
    field_idx=find(strcmpi(varargin{i},names));
    if length(field_idx)~=1
        error(['option ' varargin{i} ' unclear']);
    end
    s = setfield(s,names{field_idx},varargin{i+1});
end


