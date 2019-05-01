function Y=cat_vol_ctype(Y,type)
% ______________________________________________________________________
% Y=cat_vol_conv(Y[,type])
%
% Convert datatype with checking of min/max, nan, and rounding for 
% [u]int[8|16], single, double, and char. Default round type is 'uint8'. 
% Y has to be a matrix or cell. 
%
% This function is only writen for our private use, mostly to convert 
% single to uint8. I did not check for special behavior, for extremly 
% high values or special rounding issues, or converting to larger 
% classes etc.!
% ______________________________________________________________________
%
% Robert Dahnke 
% Structural Brain Mapping Group
% University Jena
%
% $Id: cat_vol_ctype.m 1129 2017-05-09 12:58:31Z gaser $
% ______________________________________________________________________
  
  types = {'int8','int16','int32','int64','single',...
           'uint8','uint16','uint32','uint64','double'};

  if ~exist('type','var');
    type = 'uint8';
  else
    type  = lower(type);
    % check for additional information such as -le and cut name 
    ind = strfind(type,'-');
    if ~isempty(ind)
      type = type(1:ind-1);
    end
    if all(cellfun('isempty',strfind(types,type)))
      error('MATLAB:SPM:CAT:cat_vol_ctype:UnknownType', ...
            ['ERROR: cat_vol_ctype: unknown data type ''%s'' ' ...
             '(only [u]int[8|16], single, and double).'],type);
    end
  end

  
  if iscell(Y)
    % recall function 
    for yi=1:numel(Y)
      Y{yi}=cat_vol_ctype(Y{yi},type);
    end
  else
    type = types{find(~cellfun('isempty',strfind(types,type)),1,'first')};
 
    % prepare conversion
    if ~isempty(strfind(type,'int')) || ~isempty(strfind(type,'char'))
      switch class(Y)
        case {'single','double'}
          Y = single(Y);
          Y(isnan(Y)) = 0;
          Y = round(min(single(intmax(type)),max(single(intmin(type)),Y)));
        otherwise
          % this is not working for very old matlab versions
          try
            Y = int64(Y);
            Y = round(min(int64(intmax(type)),max(int64(intmin(type)),Y)));
          catch
            Y = eval([type '(Y)']);
            Y = int64(round(min(intmax(type),max(intmin(type),Y))));
          end
      end
    elseif ~isempty(strfind(type,'single'))
      Y = min(single(realmax(type)),max(single(realmin(type)),Y));
    elseif ~isempty(strfind(type,'double'))
      Y = min(double(realmax(type)),max(double(realmin(type)),Y));
    end
    
    % convert
    eval(sprintf('Y = %s(Y);',type));
  end
  
end