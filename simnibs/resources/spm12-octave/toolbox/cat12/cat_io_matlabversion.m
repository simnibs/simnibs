function year=cat_io_matlabversion(varargin)
  vers = version;
  [t1,t2,t3,year] = regexp(vers,'\(R.....\)');
  switch year{1}(end-1)
    case 'a',  year = [year{1}(3:end-2) '1'];
    case 'b',  year = [year{1}(3:end-2) '2'];
    otherwise, year = [year{1}(3:end-2) '0'];
  end
  year = str2double(year);
  
  if nargin==1
    year = varargin{1}==year;
  end
  if nargin==2
    year = varargin{1}<=year && year<=varargin{1};
  end
end