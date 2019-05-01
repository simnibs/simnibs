function S=cat_io_mergeStruct(S,SN,ri)
% _________________________________________________________________________
% Merge to structures 'S' and 'SN'.
%
%   S = cat_io_mergeStruct(S,SN)
%   
% WARNING:
%   This function is still in developent! Be careful by using it, due to
%   unexpected behaviour. Updating of structures is a complex topic with 
%   many subcases and here only a simple alignment is used!
%
% _________________________________________________________________________
% $Id: cat_io_mergeStruct.m 966 2016-07-20 14:40:20Z gaser $

  % check input
  maxri = 20; 
  if ~exist('ri','var'), ri=0; end
  
  SFN  = fieldnames(S);
  SNFN = fieldnames(SN);

  NSFN = setdiff(SNFN,SFN);
  for ni = 1:numel(NSFN)
   if isnumeric(SN(1).(NSFN{ni}))
     S(1).(NSFN{ni}) = [];
   elseif ischar(SN(1).(NSFN{ni}))
     S(1).(NSFN{ni}) = '';
   elseif isstruct(SN(1).(NSFN{ni}))
     if ri<maxri
       Stmp = cat_io_mergeStruct(struct(),SN(1).(NSFN{ni})(1),ri+1);
     else
        Stmp = struct(); 
     end
     S(1).(NSFN{ni}) = Stmp(1);
   elseif iscell(SN(1).(NSFN{ni}))
     S(1).(NSFN{ni}) = {};
   end
  end

  NSSFN = setdiff(SFN,SNFN);
  for ni = 1:numel(NSSFN)
    if isnumeric(S(1).(NSSFN{ni}))
      SN(1).(NSSFN{ni}) = [];
    elseif ischar(S(1).(NSFN{ni}))
      SN(1).(NSSFN{ni}) = '';
    elseif isstruct(S(1).(NSSFN{ni}))
      if ri<maxri
        Stmp = cat_io_mergeStruct(struct(),S(1).(NSSFN{ni})(1),ri+1);
      else
        Stmp = struct(); 
      end
      SN(1).(NSSFN{ni}) = Stmp(1);
    elseif iscell(SN(1).(NSSFN{ni}))
      SN(1).(NSSFN{ni}) = {};
    end
  end

  S = orderfields(S);
  SN = orderfields(SN);
  
  ns = numel(S);
  for sni=1:numel(SN)
    try
      S(ns + sni) = SN(sni);
    catch
      1==1;
    end
  end
end