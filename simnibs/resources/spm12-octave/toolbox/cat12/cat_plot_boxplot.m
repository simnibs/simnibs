function [out,s] = cat_plot_boxplot(data,opt)
% _________________________________________________________________________
%
% usage: vargout = cat_plot_boxplot(data,opt);
%
%  opt.notched     = 0;             % thinner at median [0 1] with 1=0.5
%  opt.symbol      = '+o';          % outlier symbols
%  opt.vertical    = 1;             % boxplot orientation 
%  opt.maxwhisker  = 1.5;           % 
%  opt.sort        = 0;             % no sorting
%                  = 1;             % sort groups (ascending)
%                  = 2;             % sort groups (descending)[inactive]
%                  = [index];       % or by a index matrix
%  opt.names       = {};            % cell of group names
%  opt.fill        = 1;             % filling of boxes
%  opt.groupnum    = 1;             % add number of elements
% [opt.groupmin    = 5;]            % minimum number of non-nan-elements
%                                     in a group [inactive]
%  opt.ylim        = [-inf inf];    % y-axis scaling
%  opt.ygrid       = 1;             % activate y-grid-lines
%  opt.box         = 1;             % plot box
%  opt.outliers    = 1;             % plot outliers
%  opt.violin      = 0;             % violin-plot: 0 - box plot; 1 - violin plot; 2 - violin + thin box plot
%  opt.boxwidth    = 0.8;           % width of box
%  opt.groupcolor  = [R G B];       % matrix with (group)-bar-color(s) 
%                                     use jet(numel(data)) 
%                                     or other color functions
%  opt.fontsize    = [];            % axis fontsize 
%                                     important for ygrid size!
%  opt.showdata    = 0;             % show data points: 0 - no; 1 - as points; 2 - as short lines (barcode plot)
%  opt.median      = 2;             % show median: 0 - no; 1 - line; 2 - with different fill colors 
%  opt.edgecolor   = 'none';        % edge color of box 
%  opt.trans       = 0.25;          % transparency of the box
%  opt.sat         = 0.50;          % satuation of the box
%
% The box plot is a graphical display that simultaneously describes several 
% important features of a data set, such as center, spread, departure from 
% symmetry, and identification of observations that lie unusually far from
% the bulk of the data.
%
% data is a matrix with one column for each dataset, or data is a cell
% vector with one cell for each dataset.
% opt.notched = 1 produces a notched-box plot. Notches represent a robust 
% estimate of the uncertainty about the median.
% opt.notched = 0 (default) produces a rectangular box plot. 
% opt.notched in (0,1) produces a notch of the specified depth.
% opt.notched values outside [0,1] are amusing if not exactly practical.
% opt.notched sets the notched for the outlier values, default notched for
% points that lie outside 3 times the interquartile range is 'o',
% default opt.notched for points between 1.5 and 3 times the interquartile
% range is '+'. 
%
% Examples
% opt.notched = '.' points between 1.5 and 3 times the IQR is marked with
% '.' and points outside 3 times IQR with 'o'.
% opt.notched = ['x','*'] points between 1.5 and 3 times the IQR is marked with
% 'x' and points outside 3 times IQR with '*'.
% opt.vertical = 0 makes the boxes horizontal, by default opt.vertical = 1.
% maxwhisker defines the length of the wh^1iskers as a function of the IQR
% (default = 1.5). If maxwhisker = 0 then boxplot displays all data  
% values outside the box using the plotting opt.notched for points that lie
% outside 3 times the IQR.   
%
% The returned matrix s has one column for each dataset as follows:
%
%    1  minimum
%    2  1st quartile
%    3  2nd quartile (median)
%    4  3rd quartile
%    5  maximum
%    6  lower confidence limit for median
%    7  upper confidence limit for median
%
% Example
%
%   title("Grade 3 heights");
%   tics("x",1:2,["girls";"boys"]);
%   axis([0,3]);
%   cat_plot_boxplot({randn(10,1)*5+140, randn(13,1)*8+135});
%
% _________________________________________________________________________
%
% Author: Alberto Terruzzi <t-albert@libero.it>
% Version: 1.4
% Created: 6 January 2002
% Copyright (C) 2002 Alberto Terruzzi
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
%
% modified by Christian Gaser (christian.gaser@uni-jena.de) and
% Robert Dahnke (robert.dahnke@uni-jena.de)
% original version was written for octave by Alberto Terruzzi
% _________________________________________________________________________
% $Id: cat_plot_boxplot.m 1247 2017-12-13 19:48:18Z gaser $

if nargin==0, help cat_plot_boxplot; return; end

% data has to be defined as cell and shouldmbe converted if numeric
if isnumeric(data)
  if size(data,1) < size(data,2)
    data = data';
  end
  sz = size(data);
  tmp = data; clear data
  data = cell(sz(2),1);
  for i = 1:sz(2)
    data{i} = tmp(:,i);
  end
end

% default parameter
if ~exist('opt','var'), opt = struct(''); end
def.notched     = 0;
def.symbol      = '+o';
def.vertical    = 1;
def.maxwhisker  = 1.5;
def.sort        = 0; 
def.names       = num2str( (1:numel(data))' );
def.fill        = 1;
def.groupcolor  = jet(numel(data));
def.symbolcolor = 'r';
def.groupnum    = 0;
def.groupmin    = 5;
def.ylim        = [];
def.ygrid       = 0;  
def.boxwidth    = 0.8;  
def.box         = 1;
def.outliers    = 1;
def.violin      = 0;  
def.fontsize    = []; % empty = default font size
def.showdata    = 0;  
def.median      = 2;          % show median: 0 - no; 1 - line; 2 - with different fill colors 
def.edgecolor   = 'none';     % edgecolor of boxes
def.changecolor = 0;          % use brighter values for double color entries, e.g. 
                              % [red red blue blue] becomes [red light-red blue light-blue] 
def.trans       = 0.25;       % transparency of boxes
def.sat         = 0.50;       % saturation of boxes
def.subsets     = false(1,numel(data)); 
def.hflip       = 0;          % glip x-axis in case of horizontal bars

opt             = cat_io_checkinopt(opt,def);
opt.notched     = max(0,min(1,opt.notched));
opt.trans       = max(0,min(1,opt.trans * (opt.sat*4) ));
if max(opt.subsets)>1
  subsets = zeros(1,numel(data)); 
  subsets(opt.subsets+1) = 1; 
  opt.subsets = mod(cumsum(subsets),2); 
end

% first data box on the top
if ~opt.vertical 
  if opt.hflip
    for di=1:numel(data), data{di} = -data{di}; end 
    opt.ylim = fliplr(-opt.ylim); 
    if isfield(opt,'ytick'), opt.ytick = fliplr(-opt.ytick); end
  else
    data           = flipud(data); 
    opt.names      = flipud(opt.names);
    opt.groupcolor = flipud(opt.groupcolor); 
    opt.subsets    = flipud(opt.subsets); 
  end
end


% always use filling for this median plot option
if opt.median == 2, opt.fill = 1; end

% either violin or box plot
if opt.violin, opt.box = 0; end
if opt.box, opt.violin = 0; end

% figure out how many data sets we have
if iscell(data), 
  nc = length(data);
  for nci=1:nc, data{nci}=data{nci}(:); end
else
  if isvector(data), data = data(:); end
  nc = columns(data);
end
opt.names = cellstr(opt.names);
if numel(opt.names) < nc
  error('ERROR:cat_plot_boxplot:names','ERROR: Too short name list.'); 
end

% update colortable
if size(opt.groupcolor,1)==1
  %if size(opt.groupcolor,1)<nc 
  %  warning('WARNING:cat_plot_boxplot:groupcolor','WARNING: To short colortable.'); 
  %end
  opt.groupcolor = repmat(opt.groupcolor(1,:),numel(data),1);
end
% if the same color is used multiple times you may want to change it...
if opt.changecolor
  tmpcolor = opt.groupcolor(1,:);
  for ci = 2:size(opt.groupcolor,1)
    if all(opt.groupcolor(ci,:)==tmpcolor)
      opt.groupcolor(ci,:) = max(0,min(1,opt.groupcolor(ci-1,:) .* repmat(1 + 0.02 * opt.changecolor,1,3))); 
    else
      tmpcolor = opt.groupcolor(ci,:);
    end
  end
end
if numel(opt.sort)>1 && numel(opt.sort) ~= nc
  error('ERROR:cat_plot_boxplot:sort','ERROR: Too sort list.'); 
end

groupnr = cellfun(@(x) sum(~isnan(x)),data);
out.sortj = 1:length(data);
rmdata = zeros(1,nc);
% remove groups with too few elemnts
% ... require addaption for many other fields like names, color, ...
if 0 && opt.groupmin>0
  rmdata=cellfun('isempty',data) | groupnr<opt.groupmin;
  if numel(opt.sort)==numel(data)
    opt.sort(rmdata) = [];
  end
  opt.names(rmdata) = [];
  opt.groupcolor(rmdata) = [];
  data(rmdata) = [];
end
if isempty(data), 
  error('ERROR:cat_plot_boxplot:data','ERROR: Not enough (non-NaN) data (may change opt.groupmin).'); 
end
  

% add number of group elements
if opt.groupnum
  for ni=1:numel(opt.names)
    opt.names{ni}=sprintf('%s[%d]',opt.names{ni},groupnr(ni));
  end
end


% sort groups by their median value or a specific order
if isfield(opt,'sort') 
  if numel(opt.sort)==1 && opt.sort
  % sort by median
    mdata = zeros(1,numel(data));
    for i=1:numel(data), mdata(i) = cat_stat_nanmedian(data{i}(:)); end
    [mdata,sorti] = sort(mdata);
    clear mdata;
    data          = data(sorti);
    opt.names     = opt.names(sorti); 
    opt.groupcolor = opt.groupcolor(sorti,:);
    tmp = out.sortj(~rmdata);
    [tmp,out.sortj(~rmdata)] = sort(tmp(sorti));
  elseif ~opt.sort
    % noting to do, just to avoid an error
  elseif numel(opt.sort)==numel(data)
  % sort by given order
    sorti = opt.sort;
    data          = data(sorti);
    opt.names     = opt.names(sorti); 
    opt.groupcolor = opt.groupcolor(sorti,:);
    tmp = out.sortj(~rmdata);
    [tmp,out.sortj(~rmdata)] = sort(tmp(sorti));
  end
end  
out.sorti = 1:nc;
out.sortj = 1:nc;

if length(opt.symbol)==1, opt.symbol(2)=opt.symbol(1); end

if opt.notched==1, opt.notched=0.5; end
a = 1-opt.notched;  


%% compute statistics
% s will contain
%    1,5    min and max
%    2,3,4  1st, 2nd and 3rd quartile
%    6,7    lower and upper confidence intervals for median
s = zeros(7,nc);
box = zeros(1,nc);
whisker_x   = ones(2,1)*[1:nc,1:nc];
whisker_y   = zeros(2,2*nc);
outliers_x  = [];
outliers_y  = [];
outliers2_x = [];
outliers2_y = [];

% get maximum data size for violin plot
if opt.violin
  n2 = 0;
  for i=1:nc
    % Get the next data set from the array or cell array
    if iscell(data), col = data{i}(:);
    else col = data(:,i); end
    % estimate # of mesh points w.r.t. data size
    n2 = max(n2,ceil(log2(numel(col)))); 
  end
  F = zeros(2^n2,nc);
  U = zeros(2^n2,nc);
end

for i=1:nc
  % Get the next data set from the array or cell array
  if iscell(data), col = data{i}(:);
  else col = data(:,i); end
  
  % Skip missing data
  col(isnan(col)) = [];
  % Remember the data length
  nd = length(col);
  box(i) = nd;
  
  % estimate kernel density for violin plot
  if opt.violin
    [tmp, f, u] = kde(col,2^n2);
    f = (f/max(f)*opt.boxwidth*0.3)'; % width of violin plot is more narrow
    F(:,i) = f;
    U(:,i) = u;
  end
  
  if (nd > 1)
    % min,max and quartiles
    s(1:5,i) = [min(col) prctile(col,[25 50 75]) max(col)]';
    % confidence interval for the median
    est = 1.57*(s(4,i)-s(2,i))/sqrt(nd);
    s(6,i) = max([s(3,i)-est, s(2,i)]);
    s(7,i) = min([s(3,i)+est, s(4,i)]);
    % whiskers out to the last point within the desired inter-quartile range
    IQR = opt.maxwhisker*(s(4,i)-s(2,i));
    whisker_y(:,i) = [min(col(col >= s(2,i)-IQR)); s(2,i)];
    whisker_y(:,nc+i) = [max(col(col <= s(4,i)+IQR)); s(4,i)];
    % outliers beyond 1 and 2 inter-quartile ranges
    ol = (col < s(2,i)-IQR & col >= s(2,i)-2*IQR) | (col > s(4,i)+IQR & col <= s(4,i)+2*IQR);
    ol2 = col < s(2,i)-2*IQR | col > s(4,i)+2*IQR;
    outliers = col(ol);
    outliers2 = col(ol2);
    
    oll1 = (col < s(2,i)-IQR & col >= s(2,i)-2*IQR);
    olh1 = (col > s(4,i)+IQR & col <= s(4,i)+2*IQR);
    oll2 = col < s(2,i)-2*IQR;
    olh2 = col > s(4,i)+2*IQR;
    
    ind = 1:numel(col);
    out.names   = opt.names;
    out.indn.l1{i}  = ind(oll1);
    out.indn.l2{i}  = ind(oll2);
    out.indn.h1{i}  = ind(olh1);
    out.indn.h2{i}  = ind(olh2);
    out.matn.l1{i}  = oll1;
    out.matn.l2{i}  = oll2;
    out.matn.h1{i}  = olh1;
    out.matn.h2{i}  = olh2;

    if exist('sorti','var'), out.sorti   = sorti; end
    out.indo.l1{out.sortj(i)} = ind(oll1);
    out.indo.l2{out.sortj(i)} = ind(oll2);
    out.indo.h1{out.sortj(i)} = ind(olh1);
    out.indo.h2{out.sortj(i)} = ind(olh2);
    out.mato.l1{out.sorti(i)} = oll1;
    out.mato.l2{out.sorti(i)} = oll2;
    out.mato.h1{out.sorti(i)} = olh1;
    out.mato.h2{out.sorti(i)} = olh2;
    
    outliers_x = [outliers_x; i*ones(size(outliers))];
    outliers_y = [outliers_y; outliers];
    outliers2_x = [outliers2_x; i*ones(size(outliers2))];
    outliers2_y = [outliers2_y; outliers2];
  elseif (nd == 1)
    % all statistics collapse to the value of the point
    s(:,i) = col;
    % single point data sets are plotted as outliers.
    outliers_x = [outliers_x; i];
    outliers_y = [outliers_y; col];
  else
    % no statistics if no points
    s(:,i) = NaN;
  end
  
  if ~isempty(opt.ylim)
    outliers_y  = max(opt.ylim(1),min(opt.ylim(2),outliers_y)); 
    outliers2_y = max(opt.ylim(1),min(opt.ylim(2),outliers2_y)); 
  end
  
end

% Note which boxes don't have enough stats
chop = find(box <= 1);

% Draw a box around the quartiles, with width proportional to the number of
% items in the box. Draw notches if desired.
if opt.boxwidth<0
  box = repmat(abs(opt.boxwidth) * 0.4,1,numel(box));
else
  box = box*(opt.boxwidth/2/max(box));
end

quartile_x = ones(11,1)*[1:nc] + [-a;-1;-1;1;1;a;1;1;-1;-1;-a]*box;
quartile_y = s([3,7,4,4,7,3,6,2,2,6,3],:);

% thinner boxplot as option for violin plot
quartile_xthin = ones(11,1)*[1:nc] + [-1;-1;-1;1;1;1;1;1;-1;-1;-1]*0.035*ones(1,nc);

% quartiles below median for different filling colors
if ~opt.vertical && opt.hflip
  quartile_xl = ones(7,1)*[1:nc] + [a;1;1;-1;-1;-a;a]*box;
  quartile_yl = s([3,7,4,4,7,3,3],:);
else
  quartile_xl = ones(7,1)*[1:nc] + [a;1;1;-1;-1;-a;a]*box;
  quartile_yl = s([3,6,2,2,6,3,3],:);
end

% Draw a line through the median
median_x = ones(2,1)*[1:nc] + [-a;+a]*box;
median_y = s([3,3],:);

% Chop all boxes which don't have enough stats
quartile_x(:,chop)  = [];
quartile_y(:,chop)  = [];
quartile_xl(:,chop) = [];
quartile_yl(:,chop) = [];
whisker_x(:,[chop,chop+nc]) = [];
whisker_y(:,[chop,chop+nc]) = [];
median_x(:,chop) = [];
median_y(:,chop) = [];

% Add caps to the remaining whiskers
cap_x = whisker_x;
cap_x(1,:) = cap_x(1,:) - 0.1;
cap_x(2,:) = cap_x(2,:) + 0.1;
cap_y = whisker_y([1,1],:);
datac = cell2mat(data(:)); 
vp = 10^(1+round(abs(diff([min(datac(:)),max(datac(:))]))^(1/10) )); 

% use different limits for y-axis for violin plot
if opt.violin
  if isempty(opt.ylim) || isinf(opt.ylim(1)), opt.ylim(1) = min(U(:)); end
  if numel(opt.ylim)<2 || isinf(opt.ylim(2)), opt.ylim(2) = max(U(:)); end
else
  if isempty(opt.ylim) || isinf(opt.ylim(1)), opt.ylim(1) = floor((min(datac(:)) - abs(diff([min(datac(:)),max(datac(:))]))/10) * vp)/vp; end
  if numel(opt.ylim)<2 || isinf(opt.ylim(2)), opt.ylim(2) = ceil((max(datac(:))  + abs(diff([min(datac(:)),max(datac(:))]))/10) * vp)/vp; end
end

%% Do the plot
children0=get(gca,'Children');

qn = size(quartile_x,2);
  
if ~opt.vertical
%    tmp = median_x; median_x = median_y; median_y = tmp;
  tmp = quartile_x; quartile_x = quartile_y; quartile_y = tmp;
  tmp = quartile_xl; quartile_xl = quartile_yl; quartile_yl = tmp;
  tmp = cap_x; cap_x = cap_y; cap_y = tmp;
  tmp = whisker_x; whisker_x = whisker_y; whisker_y = tmp;
  tmp = outliers_x; outliers_x = outliers_y; outliers_y = tmp;
  tmp = outliers2_x; outliers2_x = outliers2_y; outliers2_y = tmp;

end

for i=1:qn
  
    % violin plot
    if opt.violin    
      if opt.vertical
        indn = max(find(U(:,i)<median_y(1,i)));
      else
        indn = max(find((U(:,i))<median_y(1,i)));
      end
      % correct thickness of median line
      if opt.vertical
        median_x(:,i) = [F(indn,i)+i;flipud(i-F(indn,i))];
      else
        median_x(:,i) = [F(indn,i)+i;flipud(1-F(indn,i))];
      end
      
      if opt.fill
        if opt.vertical
          if opt.trans, fill([F(:,i)+i;flipud(i-F(:,i))],[U(:,i);flipud(U(:,i))],[1 1 1],'FaceAlpha',1-opt.trans,'EdgeColor','none'); end % just a white box as background
          fill([F(:,i)+i;flipud(i-F(:,i))],[U(:,i);flipud(U(:,i))],opt.groupcolor(i,:),'FaceAlpha',opt.sat,'EdgeColor',opt.edgecolor);
        else
          if opt.trans, fill([U(:,i);flipud(U(:,i))],[F(:,i)+i;flipud(i-F(:,i))],[1 1 1],'FaceAlpha',1-opt.trans,'EdgeColor','none'); end
          fill([U(:,i);flipud(U(:,i))],[F(:,i)+i;flipud(i-F(:,i))],opt.groupcolor(i,:),'FaceAlpha',opt.sat,'EdgeColor',opt.edgecolor);
        end
        if i==1, hold on; end
        if opt.median == 2
          if opt.vertical
            if opt.trans, fill([F(1:indn,i)+i;flipud(i-F(1:indn,i))],[U(1:indn,i);flipud(U(1:indn,i))],[1 1 1],'FaceAlpha',1-opt.trans,'EdgeColor','none'); end
            fill([F(1:indn,i)+i;flipud(i-F(1:indn,i))],[U(1:indn,i);flipud(U(1:indn,i))],0.5*opt.groupcolor(i,:),'FaceAlpha',0.25,'EdgeColor','none');
          else
            if opt.trans, fill([U(1:indn,i);flipud(U(1:indn,i))],[F(1:indn,i)+i;flipud(i-F(1:indn,i))],[1 1 1],'FaceAlpha',1-opt.trans,'EdgeColor','none'); end
            fill([U(1:indn,i);flipud(U(1:indn,i))],[F(1:indn,i)+i;flipud(i-F(1:indn,i))],0.5*opt.groupcolor(i,:),'FaceAlpha',0.25,'EdgeColor','none');
          end
        end
      else
        if opt.vertical
          if opt.trans, plot([F(:,i)+i;flipud(i-F(:,i))],[U(:,i);flipud(U(:,i))],'Color',[1 1 1]); end
          plot([F(:,i)+i;flipud(i-F(:,i))],[U(:,i);flipud(U(:,i))],'Color',opt.groupcolor(i,:));
        else
          if opt.trans, plot([U(:,i);flipud(U(:,i))],[F(:,i)+i;flipud(i-F(:,i))],'Color',[1 1 1]); end
          plot([U(:,i);flipud(U(:,i))],[F(:,i)+i;flipud(i-F(:,i))],'Color',opt.groupcolor(i,:));
        end
        if i==1, hold on; end
        if opt.median == 2
          if opt.vertical
            if opt.trans, plot([F(1:indn,i)+i;flipud(i-F(1:indn,i))],[U(1:indn,i);flipud(U(1:indn,i))],[1 1 1],[1 1 1]); end
            plot([F(1:indn,i)+i;flipud(i-F(1:indn,i))],[U(1:indn,i);flipud(U(1:indn,i))],0.5*opt.groupcolor(i,:),0.5*opt.groupcolor(i,:));
          else
            if opt.trans, plot([U(1:indn,i);flipud(U(1:indn,i))],[F(1:indn,i)+i;flipud(i-F(1:indn,i))],[1 1 1],[1 1 1]); end
            plot([U(1:indn,i);flipud(U(1:indn,i))],[F(1:indn,i)+i;flipud(i-F(1:indn,i))],0.5*opt.groupcolor(i,:),0.5*opt.groupcolor(i,:));
          end
        end
      end
      
      % add thin box plot
      if opt.violin == 2
        if opt.fill
          fill(quartile_xthin(:,i), quartile_y(:,i),[0.5 0.5 0.5]); 
        else
          plot(quartile_xthin(:,i), quartile_y(:,i),'Color',[0.5 0.5 0.5]); 
        end
        plot(whisker_x(:,[i i+qn]), whisker_y(:,[i i+qn]),'Color',[0.5 0.5 0.5]);
      end
      
    end

    if opt.box
      if opt.fill
        if opt.trans, fill(quartile_x(:,i), quartile_y(:,i),'b-','FaceColor',[1 1 1],'FaceAlpha',1-opt.trans,'EdgeColor','none'); end
        fill(quartile_x(:,i), quartile_y(:,i),'b-','FaceColor',opt.groupcolor(i,:),'FaceAlpha',opt.sat,'EdgeColor','none');
        if i==1, hold on; end
        if opt.median == 2
          fill(quartile_xl(:,i), quartile_yl(:,i),'b-','FaceColor',0.5*opt.groupcolor(i,:),'FaceAlpha',0.25,'EdgeColor','none'); 
        end
      else
        plot(quartile_x(:,i), quartile_y(:,i),'Color',opt.groupcolor(i,:)); 
        if i==1, hold on; end
      end
      
      plot(cap_x(:,[i i+qn]), cap_y(:,[i i+qn]),'Color',[0.5 0.5 0.5]); 
      plot(whisker_x(:,[i i+qn]), whisker_y(:,[i i+qn]),'Color',[0.5 0.5 0.5]);
      if ~strcmp(opt.edgecolor,'none')
        plot(quartile_x(:,i), quartile_y(:,i),'Color',opt.edgecolor); 
      end
      
    end

    % optionally also show data either as points or short lines
    if opt.showdata == 1
      if opt.vertical
        plot(i*ones(1,length(data{i})),data{i}(:),'.','Color',0.25*opt.groupcolor(i,:));
      else
        plot(data{i}(:),i*ones(1,length(data{i})),'.','Color',0.25*opt.groupcolor(i,:));
      end
    elseif opt.showdata == 2
      if opt.vertical
        line(([-.025*ones(length(data{i}),1) .025*ones(length(data{i}),1)]+i)',([data{i}(:) data{i}(:)])','Color',0.25*opt.groupcolor(i,:));
      else
        line(([data{i}(:) data{i}(:)])',([-.025*ones(length(data{i}),1) .025*ones(length(data{i}),1)]+i)','Color',0.25*opt.groupcolor(i,:));
      end
    end

    if opt.median == 1
      if opt.groupcolor(i,1)>0.2 && opt.groupcolor(i,2)<0.5 && opt.groupcolor(i,3)<0.5
        plot(median_x(:,i), median_y(:,i),'Color',[0.5 0 0])
      else
        plot(median_x(:,i), median_y(:,i),'Color',[1 0 0])
      end
    end
    
  end
  
  if opt.symbol(1)~=' ' && opt.outliers
    plot(outliers_x,  outliers_y ,'MarkerSize',...
        max(4,min(8,80/nc)),'Marker',opt.symbol(1),'MarkerEdgeColor',opt.symbolcolor,'LineStyle','none')
  end
  if opt.symbol(2)~=' ' && opt.outliers
    plot(outliers2_x, outliers2_y,'MarkerSize',...
      max(4,min(8,80/nc)),'Marker',opt.symbol(2),'MarkerEdgeColor',opt.symbolcolor,'LineStyle','none');
  end
  
  % add labels
  linecolor = [0.8 0.8 0.8];
  if ~opt.vertical
    set(gca,'YTick',1:numel(opt.names),'YTickLabel',opt.names,'TickLength',[0 0],'ylim',[0.5 numel(opt.names)+0.5]);
    if ~isempty(opt.ylim)
      xlim(gca,opt.ylim);
    end
  else
    set(gca,'XTick',1:numel(opt.names),'XTickLabel',opt.names,'TickLength',[0 0],'xlim',[0.5 numel(opt.names)+0.5]);
    if ~isempty(opt.ylim)
      ylim(gca,opt.ylim);
    end
  end

  if ~isempty(opt.fontsize)
    set(gca,'FontSize',opt.fontsize);
  end

  % plot yticks
  if opt.ygrid 
    %%
    if opt.vertical, xytick = 'Ytick'; xylab = 'YTickLabel'; else xytick = 'Xtick'; xylab = 'XTickLabel'; end
    if isfield(opt,'ytick')
      ytick = opt.ytick; 
    else
      ytick=get(gca,xytick);
      if numel(ytick)<5, ytick=interp1(ytick,1:0.5:numel(ytick)); elseif numel(ytick)>10, ytick=ytick(1:2:end); end
    end
    set(gca,xytick,ytick); 
    if ~opt.vertical && opt.hflip
      set(gca,xylab,num2str(-ytick',sprintf('%%0.%df',-str2double(char(regexp(num2str(min(diff(ytick)),'%e'),'[+-]..','match'))) ) ) ); 
    else
      set(gca,xylab,num2str( ytick',sprintf('%%0.%df',-str2double(char(regexp(num2str(min(diff(ytick)),'%e'),'[+-]..','match'))) ) ) ); 
    end
    
    if ytick(1)<=opt.ylim(1)+eps,   ytick(1)=[];   end
    if ytick(end)>=opt.ylim(2)-eps, ytick(end)=[]; end
    if opt.vertical
      h1=plot(repmat([0;numel(opt.names)+1],1,numel(ytick)),[ytick;ytick],'Color',linecolor);
    else
      h1=plot([ytick;ytick],repmat([0;numel(opt.names)+1],1,numel(ytick)),'Color',linecolor);
    end
    uistack(h1,'bottom')

    % subgrid (just a brighter line without value)
    if opt.ygrid>1
      %%
      ytick1 = get(gca,xytick);
      ytick2 = ytick1(1):( diff([ytick1(1),ytick1(end)]) / ((numel(ytick1)-1)*opt.ygrid)):ytick1(end);
      ytick2 = setdiff(ytick2,ytick1);

      if opt.vertical
        h2=plot(repmat([0;numel(opt.names)+1],1,numel(ytick2)),[ytick2;ytick2],'Color',linecolor*0.25+0.75*ones(1,3));
      else
        h2=plot([ytick2;ytick2],repmat([0;numel(opt.names)+1],1,numel(ytick2)),'Color',linecolor*0.25+0.75*ones(1,3));
      end
      uistack(h2,'bottom')
    end
  end
  
  %% subsets
  if any(opt.subsets)
    f2 = [];
    if opt.vertical
      %%
      for i=find(opt.subsets)
        f2=[f2 fill([i-0.5 i+0.5 i+0.5 i-0.5],sort([ylim ylim]),'b-','FaceColor',[0 0 0],'FaceAlpha',0.04,'EdgeColor','none')]; %#ok<AGROW>
      end
    else
      %%
      for i=find(opt.subsets)
        f2=[f2 fill(sort([xlim xlim]),[i-0.5 i+0.5 i+0.5 i-0.5],'b-','FaceColor',[0 0 0],'FaceAlpha',0.04,'EdgeColor','none')]; %#ok<AGROW>
      end
    end
    for i=1:numel(f2), uistack(f2(i),'bottom'); end
    uistack(h1,'bottom')
    if opt.ygrid>1, uistack(h2,'bottom'); end
  end
   
%%
hold off

end

function y = prctile(x,p,dim)
%PRCTILE Percentiles of a sample.
%   Y = PRCTILE(X,P) returns percentiles of the values in X.  P is a scalar
%   or a vector of percent values.  When X is a vector, Y is the same size
%   as P, and Y(i) contains the P(i)-th percentile.  When X is a matrix,
%   the i-th row of Y contains the P(i)-th percentiles of each column of X.
%   For N-D arrays, PRCTILE operates along the first non-singleton
%   dimension.
%
%   Y = PRCTILE(X,P,DIM) calculates percentiles along dimension DIM.  The
%   DIM'th dimension of Y has length LENGTH(P).
%
%   Percentiles are specified using percentages, from 0 to 100.  For an N
%   element vector X, PRCTILE computes percentiles as follows:
%      1) The sorted values in X are taken as the 100*(0.5/N), 100*(1.5/N),
%         ..., 100*((N-0.5)/N) percentiles.
%      2) Linear interpolation is used to compute percentiles for percent
%         values between 100*(0.5/N) and 100*((N-0.5)/N)
%      3) The minimum or maximum values in X are assigned to percentiles
%         for percent values outside that range.
%
%   PRCTILE treats NaNs as missing values, and removes them.
%
%   Examples:
%      y = prctile(x,50); % the median of x
%      y = prctile(x,[2.5 25 50 75 97.5]); % a useful summary of x
%
%   See also IQR, MEDIAN, NANMEDIAN, QUANTILE.

%   Copyright 1993-2004 The MathWorks, Inc.
%   $Revision: 1247 $  $Date: 2017-12-13 20:48:18 +0100 (Wed, 13 Dec 2017) $


if ~isvector(p) || numel(p) == 0
    error('stats:prctile:BadPercents', ...
          'P must be a scalar or a non-empty vector.');
elseif any(p < 0 | p > 100) || ~isreal(p)
    error('stats:prctile:BadPercents', ...
          'P must take real values between 0 and 100');
end

% Figure out which dimension prctile will work along.
sz = size(x);
if nargin < 3 
    dim = find(sz ~= 1,1);
    if isempty(dim)
        dim = 1; 
    end
    dimArgGiven = false;
else
    % Permute the array so that the requested dimension is the first dim.
    nDimsX = ndims(x);
    perm = [dim:max(nDimsX,dim) 1:dim-1];
    x = permute(x,perm);
    % Pad with ones if dim > ndims.
    if dim > nDimsX
        sz = [sz ones(1,dim-nDimsX)];
    end
    sz = sz(perm);
    dim = 1;
    dimArgGiven = true;
end

% If X is empty, return all NaNs.
if isempty(x)
    if isequal(x,[]) && ~dimArgGiven
        y = nan(size(p),class(x));
    else
        szout = sz; szout(dim) = numel(p);
        y = nan(szout,class(x));
    end

else
    % Drop X's leading singleton dims, and combine its trailing dims.  This
    % leaves a matrix, and we can work along columns.
    nrows = sz(dim);
    ncols = prod(sz) ./ nrows;
    x = reshape(x, nrows, ncols);

    x = sort(x,1);
    nonnans = ~isnan(x);

    % If there are no NaNs, do all cols at once.
    if all(nonnans(:))
        n = sz(dim);
        if isequal(p,50) % make the median fast
            if rem(n,2) % n is odd
                y = x((n+1)/2,:);
            else        % n is even
                y = (x(n/2,:) + x(n/2+1,:))/2;
            end
        else
            q = [0 100*(0.5:(n-0.5))./n 100]';
            xx = [x(1,:); x(1:n,:); x(n,:)];
            y = zeros(numel(p), ncols, class(x));
            y(:,:) = interp1q(q,xx,p(:));
        end

    % If there are NaNs, work on each column separately.
    else
        % Get percentiles of the non-NaN values in each column.
        y = nan(numel(p), ncols, class(x));
        for j = 1:ncols
            nj = find(nonnans(:,j),1,'last');
            if nj > 0
                if isequal(p,50) % make the median fast
                    if rem(nj,2) % nj is odd
                        y(:,j) = x((nj+1)/2,j);
                    else         % nj is even
                        y(:,j) = (x(nj/2,j) + x(nj/2+1,j))/2;
                    end
                else
                    q = [0 100*(0.5:(nj-0.5))./nj 100]';
                    xx = [x(1,j); x(1:nj,j); x(nj,j)];
                    y(:,j) = interp1q(q,xx,p(:));
                end
            end
        end
    end

    % Reshape Y to conform to X's original shape and size.
    szout = sz; szout(dim) = numel(p);
    y = reshape(y,szout);
end
% undo the DIM permutation
if dimArgGiven
     y = ipermute(y,perm);  
end

% If X is a vector, the shape of Y should follow that of P, unless an
% explicit DIM arg was given.
if ~dimArgGiven && isvector(x)
    y = reshape(y,size(p)); 
end
end

function [bandwidth,density,xmesh,cdf]=kde(data,n,MIN,MAX)
% Reliable and extremely fast kernel density estimator for one-dimensional data;
%        Gaussian kernel is assumed and the bandwidth is chosen automatically;
%        Unlike many other implementations, this one is immune to problems
%        caused by multimodal densities with widely separated modes (see example). The
%        estimation does not deteriorate for multimodal densities, because we never assume
%        a parametric model for the data.
% INPUTS:
%     data    - a vector of data from which the density estimate is constructed;
%          n  - the number of mesh points used in the uniform discretization of the
%               interval [MIN, MAX]; n has to be a power of two; if n is not a power of two, then
%               n is rounded up to the next power of two, i.e., n is set to n=2^ceil(log2(n));
%               the default value of n is n=2^12;
%   MIN, MAX  - defines the interval [MIN,MAX] on which the density estimate is constructed;
%               the default values of MIN and MAX are:
%               MIN=min(data)-Range/10 and MAX=max(data)+Range/10, where Range=max(data)-min(data);
% OUTPUTS:
%   bandwidth - the optimal bandwidth (Gaussian kernel assumed);
%     density - column vector of length 'n' with the values of the density
%               estimate at the grid points;
%     xmesh   - the grid over which the density estimate is computed;
%             - If no output is requested, then the code automatically plots a graph of
%               the density estimate.
%        cdf  - column vector of length 'n' with the values of the cdf
%  Reference:
% Kernel density estimation via diffusion
% Z. I. Botev, J. F. Grotowski, and D. P. Kroese (2010)
% Annals of Statistics, Volume 38, Number 5, pages 2916-2957.

%
%  Example:
%           data=[randn(100,1);randn(100,1)*2+35 ;randn(100,1)+55];
%              kde(data,2^14,min(data)-5,max(data)+5);

data=data(:); %make data a column vector
if nargin<2 % if n is not supplied switch to the default
    n=2^14;
end
n=2^ceil(log2(n)); % round up n to the next power of 2;
if nargin<4 %define the default  interval [MIN,MAX]
    minimum=min(data); maximum=max(data);
    Range=maximum-minimum;
    MIN=minimum-Range/2; MAX=maximum+Range/2;
end
% set up the grid over which the density estimate is computed;
R=MAX-MIN; dx=R/(n-1); xmesh=MIN+[0:dx:R]; N=length(unique(data));
%bin the data uniformly using the grid defined above;
initial_data=histc(data,xmesh)/N;  initial_data=initial_data/sum(initial_data);
a=dct1d(initial_data); % discrete cosine transform of initial data
% now compute the optimal bandwidth^2 using the referenced method
I=[1:n-1]'.^2; a2=(a(2:end)/2).^2;
% use  fzero to solve the equation t=zeta*gamma^[5](t)
t_star=root(@(t)fixed_point(t,N,I,a2),N);
% smooth the discrete cosine transform of initial data using t_star
a_t=a.*exp(-[0:n-1]'.^2*pi^2*t_star/2);
% now apply the inverse discrete cosine transform
if (nargout>1)|(nargout==0)
    density=idct1d(a_t)/R;
end
% take the rescaling of the data into account
bandwidth=sqrt(t_star)*R;
density(density<0)=eps; % remove negatives due to round-off error
if nargout==0
    figure(1), plot(xmesh,density)
end
% for cdf estimation
if nargout>3
    f=2*pi^2*sum(I.*a2.*exp(-I*pi^2*t_star));
    t_cdf=(sqrt(pi)*f*N)^(-2/3);
    % now get values of cdf on grid points using IDCT and cumsum function
    a_cdf=a.*exp(-[0:n-1]'.^2*pi^2*t_cdf/2);
    cdf=cumsum(idct1d(a_cdf))*(dx/R);
    % take the rescaling into account if the bandwidth value is required
    bandwidth_cdf=sqrt(t_cdf)*R;
end

end
%################################################################
function  out=fixed_point(t,N,I,a2)
% this implements the function t-zeta*gamma^[l](t)
l=7;
f=2*pi^(2*l)*sum(I.^l.*a2.*exp(-I*pi^2*t));
for s=l-1:-1:2
    K0=prod([1:2:2*s-1])/sqrt(2*pi);  const=(1+(1/2)^(s+1/2))/3;
    time=(2*const*K0/N/f)^(2/(3+2*s));
    f=2*pi^(2*s)*sum(I.^s.*a2.*exp(-I*pi^2*time));
end
out=t-(2*N*sqrt(pi)*f)^(-2/5);
end



%##############################################################
function out = idct1d(data)

% computes the inverse discrete cosine transform
[nrows,ncols]=size(data);
% Compute weights
weights = nrows*exp(i*(0:nrows-1)*pi/(2*nrows)).';
% Compute x tilde using equation (5.93) in Jain
data = real(ifft(weights.*data));
% Re-order elements of each column according to equations (5.93) and
% (5.94) in Jain
out = zeros(nrows,1);
out(1:2:nrows) = data(1:nrows/2);
out(2:2:nrows) = data(nrows:-1:nrows/2+1);
%   Reference:
%      A. K. Jain, "Fundamentals of Digital Image
%      Processing", pp. 150-153.
end
%##############################################################

function data=dct1d(data)
% computes the discrete cosine transform of the column vector data
[nrows,ncols]= size(data);
% Compute weights to multiply DFT coefficients
weight = [1;2*(exp(-i*(1:nrows-1)*pi/(2*nrows))).'];
% Re-order the elements of the columns of x
data = [ data(1:2:end,:); data(end:-2:2,:) ];
% Multiply FFT by weights:
data= real(weight.* fft(data));
end

function t=root(f,N)
% try to find smallest root whenever there is more than one
N=50*(N<=50)+1050*(N>=1050)+N*((N<1050)&(N>50));
tol=10^-12+0.01*(N-50)/1000;
flag=0;
while flag==0
    try
        t=fzero(f,[0,tol]);
        flag=1;
    catch
        tol=min(tol*2,.1); % double search interval
    end
    if tol==.1 % if all else fails
        t=fminbnd(@(x)abs(f(x)),0,.1); flag=1;
    end
end
end





