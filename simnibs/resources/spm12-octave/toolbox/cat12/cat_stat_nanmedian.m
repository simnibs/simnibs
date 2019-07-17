function out = cat_stat_nanmedian(in, dim)
% ----------------------------------------------------------------------
% Median, not considering NaN values. Similar usage like median() or 
% MATLAB nanmedian of the statistic toolbox. Process input as double
% due to errors in large single arrays and set data class of "out" 
% to the data class of "in" at the end of the processing,
%
% out = cat_stat_nanmedian(in,dim)
%
% Example 1:
%   a = rand(4,6,3); 
%   a(rand(size(a))>0.5)=nan; 
%   av = cat_stat_nanmedian(a,3); 
%   am = nanmedian(a,3); % of the statistical toolbox ...
%   fprintf('%0.4f %0.4f\n',([av(:),am(:)])');
%
% Example 2 - special test call of example 1:
%   cat_stat_nanmedian('test')
%
% ----------------------------------------------------------------------
% Robert Dahnke 
% Structural Brain Mapping Group
% University Jena 
% ----------------------------------------------------------------------
% $Id: cat_stat_nanmedian.m 1036 2016-10-18 14:26:32Z dahnke $

  if nargin < 1
    help cat_stat_nanmedian;
    return;
  end;
  
  if ischar(in) && strcmp(in,'test')
    a = rand(4,6,3); 
    a(rand(size(a))>0.5)=nan; 
    av = cat_stat_nanmedian(a,3); 
    am = nanmedian(a,3); % of the statistical toolbox ...
    fprintf('%0.4f %0.4f\n',([av(:),am(:)])');
    out = nanmean(av(:) - am(:)); 
    return; 
  end
  
  if nargin < 2
    if size(in,1) ~= 1
      dim = 1;
    elseif size(in,2) ~= 1
      dim = 2;
    else 
      dim = 3; 
    end;
  end;
  
  sz  = size(in);

  if isempty(in), out = nan; return; end
  tp = class(in);
  in = double(in); % single failed in large arrays
   
  % reduce to 2d matrix
  pm = [dim:max(length(size(in)),dim) 1:dim-1];
  in = reshape(permute(in,pm),size(in,dim),prod(sz)/size(in,dim));

  in = sort(in,1);
  s  = size(in,1) - sum(isnan(in));
  
  % estimate median in loop
  out = zeros(size(s));
  for i = 1:length(s)
    if s(i)>0, out(i) = cat_stat_nanmean([in(floor((s(i)+1)/2),i),in(ceil((s(i)+1)/2),i)]); else out(i)=nan; end
  end
  
  % estimate median as matrix ... doesn't work :/
  %out = nan(size(s)); si=1:numel(s);
  %out(s>0) = nanmean( [ in(([floor((s(s>0)+1)/2);si(s>0)])'); in(ceil((s(s>0)+1)/2),si(s>0)) ] );
  
  % correct for permutation
  sz(dim) = 1; out = ipermute(reshape(out,sz(pm)),pm);
  
  eval(sprintf('out = %s(out);',tp));
end