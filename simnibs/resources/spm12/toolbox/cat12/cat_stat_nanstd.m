function out = cat_stat_nanstd(in, dim)
% ----------------------------------------------------------------------
% Standard deviation, not considering NaN values. Similar usage like 
% std() or MATLAB nanstd of the statistic toolbox. Process input 
% as double due to errors in large single arrays and set data class 
% of "out" to the data class of "in" at the end of the processing.
%
% out = cat_stat_nanstd(in,dim)
%
% Example 1:
%   a = rand(4,6,3); 
%   a(rand(size(a))>0.5)=nan; 
%   av = cat_stat_nanstd(a,3); 
%   am = nanstd(a,0,3); % of the statistical toolbox ...
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
% $Id: cat_stat_nanstd.m 1122 2017-04-03 17:09:48Z dahnke $

  if nargin < 1
    help cat_stat_nanstd;
    return;
  end;
  
  if ischar(in) && strcmp(in,'test')
    a = rand(4,6,3); 
    a(rand(size(a))>0.5)=nan; 
    av = cat_stat_nanstd(a,3); 
    am = nanstd(a,0,3); % of the statistical toolbox ...
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
  mn = cat_stat_nanmean(in,dim);
 
  dm = size(in); dm(setdiff(1:numel(dm),dim)) = 1;
  tmpmn = repmat(mn,dm); clear dm; 
  tmpmn(isnan(in(:))) = 0; 
    
  % estimate std
  out = (sum( (tmpin-tmpmn).^2 , dim) ./ max(1,(size(in,dim) - sum(isnan(in),dim))-1)).^0.5;
  out(isnan(mn))=nan;
  
  %eval(sprintf('out = %s(out);',tp));
end