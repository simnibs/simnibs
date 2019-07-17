function out = cat_stat_nanmean(in, dim)
% ----------------------------------------------------------------------
% Average, not considering NaN values. Similar usage like mean() or 
% MATLAB nanmean of the statistic toolbox. Process input as double
% due to errors in large single arrays and set data class of "out" 
% to the data class of "in" at the end of the processing.
%
% out = cat_stat_nanmean(in,dim)
%
% Example 1:
%   a = rand(4,6,3); 
%   a(rand(size(a))>0.5)=nan; 
%   av = cat_stat_nanmean(a,3); 
%   am = nanmean(a,3); % of the statistical toolbox ...
%   fprintf('%0.4f %0.4f\n',([av(:),am(:)])');
%
% Example 2 - special test call of example 1:
%   cat_stat_nanstd('test')
%
% ----------------------------------------------------------------------
% Robert Dahnke 
% Structural Brain Mapping Group
% University Jena 
% ----------------------------------------------------------------------
% $Id: cat_stat_nanmean.m 1122 2017-04-03 17:09:48Z dahnke $

  if nargin < 1
    help cat_stat_nanmean;
    return;
  end;
  
  if ischar(in) && strcmp(in,'test')
    a = rand(4,6,3); 
    a(rand(size(a))>0.5)=nan; 
    av = cat_stat_nanmean(a,3); 
    am = nanmean(a,3); % of the statistical toolbox ...
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
  
  if isempty(in), out = nan; return; end
  
  % estimate mean
  %tp    = class(in);
  tmpin = double(in); % single failed in large arrays
  tmpin(isnan(in(:))) = 0;
  if sum(~isnan(in),dim)==0
    out = nan(size(sum(tmpin, dim)));
  else
    out = sum(tmpin, dim) ./ max(eps,sum(~isnan(in),dim));
    out(sum(~isnan(in),dim)==0) = nan; 
  end
  
  %eval(sprintf('out = %s(out);',tp));
end