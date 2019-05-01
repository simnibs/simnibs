function varargout = cat_stat_check_cov(job)
%cat_stat_check_cov to check covariance across sample
%
% Images have to be in the same orientation with same voxel size
% and dimension (e.g. spatially registered images)
%
% Surfaces have to be same size (number of vertices).
%_______________________________________________________________________
% Christian Gaser
% $Id: cat_stat_check_cov.m 1278 2018-02-14 16:43:39Z gaser $

global alphaval fname H YpY YpYsorted  data_array data_array_diff pos ind_sorted ind_sorted_display mean_cov FS X mesh_detected ...
mn_data mx_data V Vchanged sample isxml sorted isscatter MD show_name bplot names_changed img img_alpha

% show data by fileorder
sorted = 0;

if nargin == 0
  error('No argument given.');
end

sample = [];
G      = [];
n_subjects = 0;
names_changed = 0;
alphaval = 0.5;

% check for global scaling
if isfield(job,'gSF')
  is_gSF = 1;
  gSF = job.gSF;
else
  is_gSF = 0;
end

% read filenames for each sample and indicate sample parameter
if isfield(job,'data_vol')
  mesh_detected = 0;
  n_samples = numel(job.data_vol);
  for i=1:n_samples
    
    if size(job.data_vol{i},1) == 1 % 4D data
      [pth,nam,ext] = spm_fileparts(char(job.data_vol{i}));
      % remove ",1" at the end
      job.data_vol{i} = fullfile(pth,[nam ext]);
    end
    
    V0 = spm_data_hdr_read(char(job.data_vol{i}));
    n_subjects = n_subjects + length(V0);
      
    if i==1, V = V0;
    else,    V = [V; V0]; end

    sample = [sample, i*ones(1,length(V0))];
  end
  sep = job.gap;
else
  mesh_detected = 1;
  n_samples = numel(job.data_surf);
  for i=1:n_samples
    V0 = spm_data_hdr_read(char(job.data_surf{i}));
    n_subjects = n_subjects + length(V0);
      
    if i==1, V = V0;
    else,    V = [V; V0]; end
    sample = [sample, i*ones(1,size(job.data_surf{i},1))];
  end
end
    
if ~isempty(job.c)
  for i=1:numel(job.c)
    G = [G job.c{i}];
  end
end

if isempty(char(job.data_xml))
  isxml = 0;
  QM_names = '';
  xml_files = [];
else
  xml_files = char(job.data_xml);
end

if ~isempty(xml_files)

  isxml = 1;
  if size(xml_files,1) ~= n_subjects
    error('XML-files must have the same number as sample size');
  end
  
  QM = ones(n_subjects,3);
  QM_names = char('Noise','Bias','Weighted overall image quality');

  spm_progress_bar('Init',n_subjects,'Load xml-files','subjects completed')
  for i=1:n_subjects
    % get basename for xml- and data files
    [pth, xml_name] = fileparts(deblank(xml_files(i,:)));
    [pth, data_name] = fileparts(V(i).fname);
    
    % remove leading 'cat_'
    xml_name = xml_name(5:end);
    
    % check for filenames
    if isempty(strfind(data_name,xml_name))
      warning('Please check file names because of deviating subject names\n: %s vs. %s\n',V(i).fname,xml_files(i,:));
    end
    
    xml = cat_io_xml(deblank(xml_files(i,:)));
    try
      QM(i,:) = [xml.qualityratings.NCR xml.qualityratings.ICR xml.qualityratings.IQR];
    catch % also try to use old version
      QM(i,:) = [xml.QAM.QM.NCR xml.QAM.QM.ICR xml.QAM.QM.rms];
    end
    spm_progress_bar('Set',i);  
  end
  spm_progress_bar('Clear');
  
end

[pth,nam] = spm_fileparts(V(1).fname);

if mesh_detected
  % load surface texture data
  Y = spm_data_read(V)';
  
  % optional global scaling
  if is_gSF
    for i=1:numel(V)
      Y(:,2) = Y(:,2)*gSF(i);
    end
  end
  
  Y(isnan(Y)) = 0;
  
else
  % voxelsize and origin
  vx =  sqrt(sum(V(1).mat(1:3,1:3).^2));
  Orig = V(1).mat\[0 0 0 1]';

  if length(V)>1 && any(any(diff(cat(1,V.dim),1,1),1))
    error('images don''t all have same dimensions')
  end
  if max(max(max(abs(diff(cat(3,V.mat),1,3))))) > 1e-8
    error('images don''t all have same orientation & voxel size')
  end
end

% add constant to nuisance parameter
if ~isempty(G)
  if size(G,1) ~= n_subjects
    G = G';
  end
  G = [ones(n_subjects,1) G];
end

% positions & font size
ws = spm('Winsize','Graphics');
FS = spm('FontSizes');

pos = struct(...
    'fig',    [10  10  1.3*ws(3) 1.1*ws(3)],... % figure
    'cbar',   [0.240 0.950 0.300 0.020],... % colorbar for correlation matrix
    'corr',   [-0.02 0.050 0.825 0.825],... % correlation matrix
    'scat',   [0.050 0.050 0.700 0.825],... % scatter plot
    'close',  [0.775 0.925 0.200 0.050],... % close button
    'show',   [0.775 0.875 0.200 0.050],... % button to show worst cases
    'boxp',   [0.775 0.820 0.200 0.050],... % button to display boxplot
    'sort',   [0.775 0.775 0.200 0.050],... % button to enable ordered matrix
    'chbox',  [0.775 0.750 0.200 0.050],... % show filenames?
    'text',   [0.775 0.600 0.200 0.150],... % textbox
    'aslider',[0.775 0.555 0.200 0.040],... % slider for alpha overlay
    'slice',  [0.775 0.050 0.200 0.400],... % two single images according to position of mouse pointer
    'sslider',[0.775 0.010 0.200 0.040]);   % slider for z-slice   

if mesh_detected
  % rescue unscaled data min/max
  mn_data = min(Y(:));
  mx_data = max(Y(:));
  Y = Y - repmat(mean(Y,2), [1 size(Y,2)]);

  % remove nuisance and add mean again (otherwise correlations are quite small and misleading)
  if ~isempty(G) 
    Ymean = repmat(mean(Y), [n_subjects 1]);
    Y = Y - G*(pinv(G)*Y) + Ymean;
  end

  data_array = Y';
  YpY = (Y*Y')/n_subjects;

  % calculate residual mean square of mean adjusted Y
  Y = Y - repmat(mean(Y,1), [n_subjects 1]);
  data_array_diff = Y';
  
  MSE = sum(Y.*Y,2);

  clear Y
else
  % consider image aspect ratio
  pos.slice = [0.775 0.050 0.20 0.40*V(1).dim(2)/V(1).dim(1)];

  slices = 1:sep:V(1).dim(3);

  dimx = length(1:sep:V(1).dim(1));
  dimy = length(1:sep:V(1).dim(2));
  Y = zeros(n_subjects, prod(dimx*dimy));
  YpY = zeros(n_subjects);
  MSE = zeros(n_subjects,1);
  data_array = zeros([V(1).dim(1:2) n_subjects]);

  %-Start progress plot
  %-----------------------------------------------------------------------
  spm_progress_bar('Init',V(1).dim(3),'Check correlation','planes completed')

  for j=slices

    M  = spm_matrix([0 0 j 0 0 0 sep sep sep]);

    for i = 1:n_subjects
      img = spm_slice_vol(V(i),M,[dimx dimy],[1 0]);
      img(isnan(img)) = 0;
      Y(i,:) = img(:);
      if is_gSF
        Y(i,:) = Y(i,:)*gSF(i);
      end
    end

    % make sure data is zero mean
    Y = Y - repmat(mean(Y,2), [1 prod(dimx*dimy)]);

    % remove nuisance and add mean again (otherwise correlations are quite small and misleading)
    if ~isempty(G) 
      Ymean = repmat(mean(Y), [n_subjects 1]);
      Y = Y - G*(pinv(G)*Y) + Ymean;
    end
    YpY = YpY + (Y*Y')/n_subjects;

    % calculate residual mean square of mean adjusted Y
    Y = Y - repmat(mean(Y,1), [n_subjects 1]);
    
    MSE = MSE + sum(Y.*Y,2);

    spm_progress_bar('Set',j);  

  end

  % correct filenames for 4D data
  if strcmp(V(1).fname, V(2).fname)
    names_changed = 1;
    Vchanged = V;
    for i=1:n_subjects
      [pth,nam,ext] = spm_fileparts(V(i).fname);
      V(i).fname = fullfile(pth, [nam sprintf('%04d',i) ext]);
    end
  end
  
  clear Y
  spm_progress_bar('Clear');
end

% normalize YpY
d      = sqrt(diag(YpY)); % sqrt first to avoid under/overflow
dd     = d*d';
YpY    = YpY./(dd+eps);
t      = find(abs(YpY) > 1); 
YpY(t) = YpY(t)./abs(YpY(t));
YpY(1:n_subjects+1:end) = sign(diag(YpY));

% extract mean correlation for each data set
mean_cov = zeros(n_subjects,1);
for i=1:n_subjects
  % extract row for each subject
  cov0 = YpY(i,:);

  % remove cov with its own
  cov0(i) = [];
  mean_cov(i) = mean(cov0);
end

fprintf('\n');
fname_m = [];
fname_tmp = cell(n_samples,1);
fname_s   = cell(n_samples,1);
fname_e   = cell(n_samples,1);

for i=1:n_samples
  [tmp, fname_tmp{i}] = spm_str_manip(char(V(sample == i).fname),'C');
  fname_m = [fname_m; fname_tmp{i}.m];
  fname_s{i} = fname_tmp{i}.s;
  fprintf('Compressed filenames sample %d: %s  \n',i,tmp);
end

fname = struct('s',{fname_s},'e',{fname_e},'m',{fname_m});

% print suspecious files with cov>0.925
YpY_tmp = YpY - tril(YpY);
[indx, indy] = find(YpY_tmp>0.925);

% if more than 25% of the data this points to longitudinal data of one subject
% and no warning will appear
if ~isempty(indx) && (sqrt(length(indx)) < 0.25*n_subjects)
  fprintf('\nUnusual large correlation (check that subjects are not identical):\n');
  for i=1:length(indx)
    % exclude diagonal
    if indx(i) ~= indy(i)
      % report file with lower mean correlation first
      if mean_cov(indx(i)) < mean_cov(indy(i))
        fprintf('%s and %s: %3.3f\n',fname.m{indx(i)},fname.m{indy(i)},YpY(indx(i),indy(i)));
      else
        fprintf('%s and %s: %3.3f\n',fname.m{indy(i)},fname.m{indx(i)},YpY(indy(i),indx(i)));
      end
    end
  end
end

% sort data
[mean_cov_sorted, ind_sorted] = sort(mean_cov,'descend');
YpYsorted = YpY(ind_sorted,ind_sorted);

ind_sorted_display = ind_sorted;

threshold_cov = mean(mean_cov) - 2*std(mean_cov);
n_thresholded = min(find(mean_cov_sorted < threshold_cov));

if ~isempty(n_thresholded)
  fprintf('\nThese data have a mean correlation below 2 standard deviations.\n');
  fprintf('This does not necessarily mean that you have to exclude these data. However, these data have to be carefully checked:\n');
  for i=n_thresholded:n_subjects
    fprintf('%s: %3.3f\n',V(ind_sorted(i)).fname,mean_cov_sorted(i));
  end
end

if nargout>0
  varargout{1} = struct('table',[cellstr(V.fname),num2cell(mean_cov)],...
                        'covmat',YpY,...
                        'sorttable',[cellstr(V(ind_sorted).fname),num2cell(mean_cov_sorted)],...
                        'sortcovmat',YpYsorted, ...
                        'cov',mean_cov,...
                        'threshold_cov',threshold_cov);
end

% create figure
H.figure = figure(2);
clf(H.figure);

set(H.figure,'MenuBar','none','Position',pos.fig,'NumberTitle','off');
    
if mesh_detected
  set(H.figure,'Name','Click in image to display surfaces');
else
  set(H.figure,'Name','Click in image to display slices');
end

cm = datacursormode(H.figure);
set(cm,'UpdateFcn',@myupdatefcn,'SnapToDataVertex','on','Enable','on');
try set(cm,'NewDataCursorOnClick',false); end

% add colorbar
H.cbar = axes('Position',pos.cbar,'Parent',H.figure);
image((1:64));

isscatter = 0;
show_matrix(YpY, sorted);

% create two colormaps
cmap = [hot(64); gray(64)];
colormap(cmap)

% display YTick with 5 values (limit accuracy for floating numbers)
set(H.cbar,'YTickLabel','','XTickLabel','','XTick',linspace(1,64,5), 'XTickLabel',...
  round(100*linspace(min(YpY(:)),max(YpY(:)),5))/100,'TickLength',[0 0]);

% add button for closing all windows
H.close = uicontrol(H.figure,...
        'string','Close','Units','normalized',...
        'position',pos.close,...
        'Style','Pushbutton','HorizontalAlignment','center',...
        'callback','for i=2:26, try close(i); end; end;',...
        'ToolTipString','Close windows',...
        'Interruptible','on','Enable','on',...
        'FontSize',FS(7));

% check button
H.show = uicontrol(H.figure,...
        'string','Check most deviating data','Units','normalized',...
        'position',pos.show,...
        'Style','Pushbutton','HorizontalAlignment','center',...
        'callback',@check_worst_data,...
        'ToolTipString','Display most deviating files',...
        'Interruptible','on','Enable','on',...
        'FontSize',FS(7));

show_name = 0;

% create popoup menu 
if isxml

  % estimate Mahalanobis distance between mean corr. and weighted overall quality
  X = [mean_cov, QM(:,3)]; % mean correlation and IQR
  S = cov(X);
  mu = mean(X);
  MD = (X-repmat(mu,[length(X),1]))*inv(S)*(X-repmat(mu,[length(X),1]))';
  MD = diag(MD);
  
  str  = { 'Boxplot...','Mean correlation',QM_names,'Mahalanobis distance'};
  tmp  = { {@show_mean_boxplot, mean_cov, 'Mean correlation  ', 1},...
           {@show_mean_boxplot, QM(:,1), QM_names(1,:), -1},...
           {@show_mean_boxplot, QM(:,2), QM_names(2,:), -1},...
           {@show_mean_boxplot, QM(:,3), QM_names(3,:), -1},...
           {@show_mean_boxplot, MD, 'Mahalanobis distance  ', -1} };
else
  str  = { 'Boxplot...','Mean correlation'};
  tmp  = { {@show_mean_boxplot, mean_cov, 'Mean correlation  ', 1} };
end

H.boxp = uicontrol(H.figure,...
        'string',str,'Units','normalized',...
        'position',pos.boxp,'UserData',tmp,...
        'Style','PopUp','HorizontalAlignment','center',...
        'callback','spm(''PopUpCB'',gcbo)',...
        'ToolTipString','Display boxplot',...
        'Interruptible','on','Visible','on',...
        'FontSize',FS(7));

if isxml
  str  = { 'Image...','Mean Correlation: Order by selected filenames','Mean Correlation: Sorted by mean correlation','Mahalanobis distance'};
  tmp  = { {@show_matrix, YpY, 0},...
           {@show_matrix, YpYsorted, 1},...
           {@show_mahalanobis, X} };
else
  str  = { 'Correlation matrix...','Order by selected filename','Sorted by mean correlation'};
  tmp  = { {@show_matrix, YpY, 0},...
           {@show_matrix, YpYsorted, 1} };
end

H.sort = uicontrol(H.figure,...
        'string',str,'Units','normalized',...
        'position',pos.sort,'UserData',tmp,...
        'Style','PopUp','HorizontalAlignment','center',...
        'callback','spm(''PopUpCB'',gcbo)',...
        'ToolTipString','Sort matrix',...
        'Interruptible','on','Visible','on',...
        'FontSize',FS(7));

H.chbox = uicontrol(H.figure,...
        'string','Show filenames in boxplot','Units','normalized',...
        'position',pos.chbox,...
        'Style','CheckBox','HorizontalAlignment','center',...
        'callback',{@checkbox_names},...
        'ToolTipString','Show filenames in boxplot',...
        'Interruptible','on','Visible','on',...
        'BackgroundColor',[0.8 0.8 0.8],'FontSize',FS(6));

H.text = uicontrol(H.figure,...
        'Units','normalized','position',pos.text,...
        'String','Click in image to display slices',...
        'Style','text','HorizontalAlignment','center',...
        'ToolTipString','Select slice for display',...
        'FontSize',FS(6));

H.alpha = uicontrol(H.figure,...
        'Units','normalized','position',pos.aslider,...
        'Min',0,'Max',1,...
        'Style','slider','HorizontalAlignment','center',...
        'callback',@update_alpha,'Value',0.5,...
        'ToolTipString','Change Opacity of pos. (green colors) and neg. (red colors) differences to sample mean',...
        'SliderStep',[0.01 0.1],'Visible','off');

H.alpha_txt = uicontrol(H.figure,...
        'Units','normalized','HorizontalAlignment','center',...
        'Style','text','BackgroundColor',[0.8 0.8 0.8],...
        'Position',[pos.aslider(1) pos.aslider(2)-0.005 0.2 0.02],...
        'String','Overlay Opacity of Differences to Sample Mean',...
        'FontSize',FS(6),'Visible','off');

% add slider only for volume data
if ~mesh_detected
  H.mm = uicontrol(H.figure,...
        'Units','normalized','position',pos.sslider,...
        'Min',(1 - Orig(3))*vx(3),'Max',(V(1).dim(3) - Orig(3))*vx(3),...
        'Style','slider','HorizontalAlignment','center',...
        'callback',@update_slices_array,...
        'ToolTipString','Select slice for display',...
        'SliderStep',[0.005 0.05],'Visible','off');

  H.mm_txt = uicontrol(H.figure,...
        'Units','normalized','HorizontalAlignment','center',...
        'Style','text','BackgroundColor',[0.8 0.8 0.8],...
        'Position',[pos.sslider(1) pos.sslider(2)-0.005 0.2 0.02],...
        'String','Slice [mm]','Visible','off','FontSize',FS(6));
        
  update_slices_array;
end
   
show_mean_boxplot(mean_cov,'Mean correlation  ',1);

% check for replicates
for i=1:n_subjects
  for j=1:n_subjects
    if (i>j) && (mean_cov(i) == mean_cov(j))
      try
        nami = deblank(V(i).fname);
        namj = deblank(V(j).fname);
        if strcmp(nami(end-1:end),',1')
          nami = nami(1:end-2);
        end 
        if strcmp(namj(end-1:end),',1')
          namj = namj(1:end-2);
        end 
        s = unix(['diff ' nami ' ' namj]);
        if (s==0), fprintf(['\nWarning: ' nami ' and ' namj ' are same files?\n']); end
      end
    end
  end
end

%-End
%-----------------------------------------------------------------------

%-----------------------------------------------------------------------
function check_worst_data(obj, event_obj)
%-----------------------------------------------------------------------
global V ind_sorted_display mesh_detected mn_data mx_data H

n = length(V);
number = min([n 24]);
number = spm_input('How many files ?',1,'e',number);
number = min([number 24]);
number = min([number length(V)]);
  
list = char(V(ind_sorted_display(n:-1:1)).fname);
list2 = list(1:number,:);

if mesh_detected
  % display single meshes and correct colorscale of colorbar
  for i=1:number
    h = cat_surf_render('Disp',deblank(list2(i,:)));
    
    % shift each figure slightly
    if i==1
        pos = get(h.figure,'Position');
    else
        pos = pos - [20 20 0 0];
    end
    
    % remove menubar and toolbar, use filename as title
    set(h.figure,'MenuBar','none','Toolbar','none','Name',spm_file(list2(i,:),'short50'),...
         'NumberTitle','off','Position',pos);
    cat_surf_render('ColourMap',h.axis,jet);
    cat_surf_render('ColourBar',h.axis,'on');
    cat_surf_render('CLim',h,[mn_data mx_data]);
  end
else
  spm_check_registration(list2);
  spm_orthviews('Resolution',0.2);
  set(H.boxp,'Visible','on');
end
return

%-----------------------------------------------------------------------
function checkbox_names(obj, event_obj)
%-----------------------------------------------------------------------
global H show_name data_boxp name_boxp quality_order

  show_name = get(H.chbox,'Value');
  show_mean_boxplot;
  
return
        
%-----------------------------------------------------------------------
function show_mahalanobis(X)
%-----------------------------------------------------------------------
global H FS pos isscatter ind_sorted_display MD

% clear larger area and set background color to update labels and title
H.ax = axes('Position',[-.1 -.1 1.1 1.1],'Parent',H.figure);
cla(H.ax);
set(H.ax,'Color',[0.8 0.8 0.8]);

H.ax = axes('Position',pos.scat,'Parent',H.figure);

cmap = [jet(64); gray(64)];

S = cov(X);
mu = mean(X);
MD = (X-repmat(mu,[length(X),1]))*inv(S)*(X-repmat(mu,[length(X),1]))';
MD = diag(MD);

% because we use a splitted colormap we have to set the color
% values explicitely
MD2 = 63*MD/max(MD);
C = zeros(length(MD),3);
for i=1:length(MD)
  C(i,:) = cmap(round(MD2(i))+1,:);
end
scatter(X(:,1),X(:,2),30,C,'o','Linewidth',2);

xlabel('<----- Worst ---      Mean correlation      --- Best ------>  ','FontSize',FS(8),'FontWeight','Bold');
ylabel('<----- Best ---      Weighted overall image quality      --- Worst ------>  ','FontSize',FS(8),'FontWeight','Bold');
title('<--- Best -- Mahalanobis distance -- Worst ---->  ','FontSize',FS(10),'FontWeight','Bold');

% add colorbar
H.cbar = axes('Position',pos.cbar,'Parent',H.figure);
image((1:64));

% display YTick with 5 values (limit accuracy for floating numbers)
set(H.cbar,'YTickLabel','','XTickLabel','','XTick',linspace(1,64,5), 'XTickLabel',...
  round(100*linspace(min(MD),max(MD),5))/100,'TickLength',[0 0]);

% update index of worst files
[tmp, ind_sorted_display] = sort(MD,'ascend');

colormap(cmap)

isscatter = 1;

return

%-----------------------------------------------------------------------
function show_matrix(data, order)
%-----------------------------------------------------------------------
global H FS pos sorted YpY isscatter mesh_detected

% get sorting order
sorted = order;

% clear larger area and set background color to update labels and title
H.ax = axes('Position',[-.1 -.1 1.1 1.1],'Parent',H.figure);
cla(H.ax);
set(H.ax,'Color',[0.8 0.8 0.8]);

H.ax = axes('Position',pos.corr,'Parent',H.figure);

% scale data to 0..1
mn = min(data(:));
mx = max(data(:));
data_scaled = (data - mn)/(mx - mn);

% show only lower left triangle
ind_tril = find(tril(ones(size(data))));
ima = zeros(size(data));
ima(ind_tril) = data_scaled(ind_tril);
image(64*ima)
set(gca,'XTickLabel','','YTickLabel','');
axis image

if sorted
  xlabel('<----- Best ---      File Order      --- Worst ------>  ','FontSize',FS(8),'FontWeight','Bold');
  ylabel('<----- Worst ---      File Order      --- Best ------>  ','FontSize',FS(8),'FontWeight','Bold');
  title('Sorted Sample Correlation Matrix  ','FontSize',FS(10),'FontWeight','Bold');
else
  xlabel('<----- First ---      File Order      --- Last ------>  ','FontSize',FS(8),'FontWeight','Bold');
  ylabel('<----- Last ---      File Order      --- First ------>  ','FontSize',FS(8),'FontWeight','Bold');
  title('Sample Correlation Matrix  ','FontSize',FS(10),'FontWeight','Bold');
end

H.cbar = axes('Position',pos.cbar,'Parent',H.figure);
image((1:64));

% display YTick with 5 values (limit accuracy for floating numbers)
set(H.cbar,'YTickLabel','','XTickLabel','','XTick',linspace(1,64,5), 'XTickLabel',...
  round(100*linspace(min(YpY(:)),max(YpY(:)),5))/100,'TickLength',[0 0]);

cmap = [hot(64); gray(64)];
colormap(cmap)

isscatter = 0;

% remove sliders and text
try
  if ~mesh_detected
    set(H.mm,'Visible','off'); 
    set(H.mm_txt,'Visible','off');
  end
  set(H.alpha,'Visible','off');
  set(H.alpha_txt,'Visible','off');
end

return

%-----------------------------------------------------------------------
function show_mean_boxplot(data_boxp, name_boxp, quality_order)
%-----------------------------------------------------------------------
global fname FS sample ind_sorted_display show_name bp

if nargin == 0
  data_boxp = bp.data;
  name_boxp = bp.name;
  quality_order = bp.order;
end

Fgraph = spm_figure('GetWin','Graphics');
spm_figure('Clear',Fgraph);
set(Fgraph,'Renderer','OpenGL');

n_samples = max(sample);

xpos = cell(1,n_samples);
data = cell(1,n_samples);

allow_violin = 2;

hold on
for i=1:n_samples
  ind = find(sample == i);
  data{i} = data_boxp(ind);
  
  if length(ind) < 8
    allow_violin = 0;
  end
  
  if n_samples == 1
    xpos{i} = (i-1)+2*(0:length(ind)-1)/(length(ind)-1);
  else
    xpos{i} = 0.5/length(ind) + 0.5+(i-1)+1*(0:length(ind)-1)/(length(ind));
  end

  for j=1:length(ind)
    if show_name
      text(xpos{i}(j),data{i}(j),fname.m{ind(j)},'FontSize',FS(7),'HorizontalAlignment','center')
    else
      plot(xpos{i}(j),data{i}(j),'k.');
    end
  end
end

opt = struct('groupnum',0,'ygrid',0,'violin',allow_violin,'median',2,'groupcolor',jet(n_samples));
ylim_add = 0.075;

cat_plot_boxplot(data,opt);

set(gca,'XTick',[],'XLim',[-.25 n_samples+1.25]);
if max(data_boxp) > min(data_boxp)
  yamp = max(data_boxp) - min(data_boxp);
  ylim_min = min(data_boxp) - ylim_add*yamp;
  ylim_max = max(data_boxp) + ylim_add*yamp;
  set(gca,'YLim',[ylim_min ylim_max]);
end

% add colored labels and title
if n_samples > 1
  [tmp,  tmp2] = spm_str_manip(char(fname.s),'C');
  title_str = sprintf('Boxplot: %s  \n%s ',name_boxp, strrep(tmp,tmp2.s,''));
  fprintf('\nCommon filename: %s\n',tmp);
else
  title_str = sprintf('Boxplot: %s  \nCommon filename: %s*',name_boxp,spm_file(char(fname.s),'short25'));
end
title(title_str,'FontSize',FS(8),'FontWeight','Bold');
xlabel('<----- First ---      File Order      --- Last ------>  ','FontSize',FS(10),...
    'FontWeight','Bold');

xpos = -0.40 - n_samples*0.1;

if (length(data_boxp) > 2)
  if quality_order > 0
    text(xpos, ylim_min,'<----- Low rating (poor quality)  ','Color','red','Rotation',...
        90,'HorizontalAlignment','left','FontSize',FS(9),'FontWeight','Bold')
    text(xpos, ylim_max,'High rating (good quality) ------>  ','Color','green','Rotation',...
        90,'HorizontalAlignment','right','FontSize',FS(9),'FontWeight','Bold')
  else
    text(xpos, ylim_max,'Low rating (poor quality) ------>  ','Color','red','Rotation',...
        90,'HorizontalAlignment','right','FontSize',FS(9),'FontWeight','Bold')
    text(xpos, ylim_min,'<----- High rating (good quality)  ','Color','green','Rotation',...
        90,'HorizontalAlignment','left','FontSize',FS(9),'FontWeight','Bold')
  end
  text(xpos, (ylim_max+ylim_min)/2,sprintf('%s',name_boxp),'Color','black','Rotation',...
        90,'HorizontalAlignment','center','FontSize',FS(9),'FontWeight','Bold')
end

hold off

% estimate sorted index new for displaying worst files
if quality_order > 0
  [tmp, ind_sorted_display] = sort(data_boxp,'descend');
else
  [tmp, ind_sorted_display] = sort(data_boxp,'ascend');
  
end

bp = struct('data',data_boxp,'name',name_boxp,'order',quality_order);

return

%-----------------------------------------------------------------------
function update_alpha(obj, event_obj)
%-----------------------------------------------------------------------
global alphaval H mesh_detected img img_alpha

if isfield(H,'alpha')
  alphaval = get(H.alpha,'Value');
else
  alphaval = 0.5;
end

% display image with 2nd colorbar (gray)
image(65 + img);
if ~mesh_detected, axis image; end
set(gca,'XTickLabel','','YTickLabel','','TickLength',[0 0]);

% prepare alpha overlays for red and green colors
if alphaval > 0
  % get 2%/98% ranges of difference image
  range = cat_vol_iscaling(img_alpha(:),[0.02 0.98]);

  hold on
  alpha_g = cat(3, zeros(size(img_alpha)), alphaval*ones(size(img_alpha)), zeros(size(img_alpha)));
  alpha_r = cat(3, alphaval*ones(size(img_alpha)), zeros(size(img_alpha)), zeros(size(img_alpha)));
  hg = image(alpha_g); set(hg, 'AlphaData', img_alpha.*(img_alpha>range(2)),'AlphaDataMapping','scaled')
  if ~mesh_detected, axis image; end
  hr = image(alpha_r); set(hr, 'AlphaData',-img_alpha.*(img_alpha<range(1)),'AlphaDataMapping','scaled')
  if ~mesh_detected, axis image; end
  hold off
end

return

%-----------------------------------------------------------------------
function update_slices_array(obj, event_obj)
%-----------------------------------------------------------------------
global alphaval V Vchanged fname data_array data_array_diff H YpY pos sorted ind_sorted isscatter names_changed FS img_alpha img

if isfield(H,'mm')
  slice_mm = get(H.mm,'Value');
else
  slice_mm = 0;
end

if names_changed
  P = Vchanged;
else
  P = V;
end

vx   =  sqrt(sum(P(1).mat(1:3,1:3).^2));
Orig = P(1).mat\[0 0 0 1]';
sl   = round(slice_mm/vx(3)+Orig(3));

% if slice is outside of image use middle slice
if (sl>P(1).dim(3)) || (sl<1)
  sl = round(P(1).dim(3)/2);
end

M  = spm_matrix([0 0 sl]);
data_array_diff = data_array;

for i = 1:length(V)
  img = spm_slice_vol(P(i),M,P(1).dim(1:2),[1 0]);
  img(isnan(img)) = 0;
  
  % rescue unscaled data
  data_array_diff(:,:,i) = img;

  % scale image according to mean
  data_array(:,:,i) = img/mean(img(img ~= 0));
end

% calculate individual difference to mean image
for i=1:size(data_array_diff,3)
  data_array_diff(:,:,i) = data_array_diff(:,:,i) - mean(data_array_diff,3);
end

% enhance contrast and scale image to 0..64
mn = min(data_array(:));
mx = max(data_array(:));
data_array = 64*((data_array - mn)/(mx-mn));

if sorted
  if isfield(pos,'x')
    x = ind_sorted(pos.x);
    if ~isscatter
      y = ind_sorted(pos.y);
    end
  end
else
  if isfield(pos,'x')
    x = pos.x;
    if ~isscatter
      y = pos.y;
    end
  end
end

% check whether mouse position is defined
if isfield(pos,'x')
  if isscatter
    img       = data_array(:,:,x)';
    img_alpha = data_array_diff(:,:,x)';
  else
    img       = [data_array(:,:,y) data_array(:,:,x)]';
    img_alpha = [data_array_diff(:,:,y) data_array_diff(:,:,x)]';
  end
  
  % correct orientation
  img = rot90(img,2);
  img_alpha = rot90(img_alpha,2);
  
  % use gray scale colormap for values > 64
  axes('Position',pos.slice);
  image(65 + img);
  axis image
  set(gca,'XTickLabel','','YTickLabel','');
  
  % prepare alpha overlays for red and green colors
  if alphaval > 0
    % get 2%/98% ranges of difference image
    range = cat_vol_iscaling(img_alpha(:),[0.02 0.98]);

    hold on
    alpha_g = cat(3, zeros(size(img_alpha)), alphaval*ones(size(img_alpha)), zeros(size(img_alpha)));
    alpha_r = cat(3, alphaval*ones(size(img_alpha)), zeros(size(img_alpha)), zeros(size(img_alpha)));
    hg = image(alpha_g); set(hg, 'AlphaData', img_alpha.*(img_alpha>range(2)),'AlphaDataMapping','scaled')
    axis image
    hr = image(alpha_r); set(hr, 'AlphaData',-img_alpha.*(img_alpha<range(1)),'AlphaDataMapping','scaled')
    axis image
    hold off
  end
  
  if isscatter
    txt = {sprintf('%s',spm_file(fname.m{x},'short25')),[],['Displayed slice: ',num2str(round(get(H.mm,'Value'))),' mm']};
  else
    txt = {sprintf('Correlation: %3.3f',YpY(x,y)),[],['Top: ',spm_file(fname.m{x},'short25')],...
      ['Bottom: ',spm_file(fname.m{y},'short25')],[],['Displayed slice: ',num2str(round(get(H.mm,'Value'))),' mm']};
  end
  set(H.text,'String',txt,'FontSize',FS(6));
  set(H.mm_txt,'String',[num2str(round(get(H.mm,'Value'))),' mm'],...
      'FontSize',FS(6));
end

return

%-----------------------------------------------------------------------
function txt = myupdatefcn(obj, event_obj)
%-----------------------------------------------------------------------
global alphaval fname sample H X YpY data_array data_array_diff pos mesh_detected ind_sorted sorted isscatter img img_alpha FS

set(H.alpha,'Visible','on');
set(H.alpha_txt,'Visible','on');

pos_mouse = get(event_obj, 'Position');

if isscatter
  pos.x = find(X(:,1) == pos_mouse(1));
  if isempty(pos.x)
    pos.x = find(X(:,2) == pos_mouse(2));
  end

  % text info for data cursor window
  txt = {sprintf('%s',fname.m{pos.x})};

  % text info for textbox
  txt2 = {sprintf('%s',spm_file(fname.m{pos.x},'short25')),[],'Difference to Sample Mean (red: - green: +)'};

  set(H.text,'String',txt2,'FontSize',FS(6));
  axes('Position',pos.slice);

  x = pos.x;
else
  % check for valid mouse position
  if pos_mouse(1) > pos_mouse(2) || pos_mouse(1)>length(sample) || pos_mouse(2)>length(sample)
    txt = {''};
    return
  end

  % save position of mouse
  pos.x = pos_mouse(1);
  pos.y = pos_mouse(2);

  if sorted
    if isfield(pos,'x')
      x = ind_sorted(pos.x);
      y = ind_sorted(pos.y);
    end
  else
    if isfield(pos,'x')
      x = pos.x;
      y = pos.y;
    end
  end

  % text info for data cursor window
  if mesh_detected
    txt = {sprintf('Correlation: %3.3f',YpY(x,y)),['Left: ',fname.m{x}],...
      ['Right: ',fname.m{y}]};
  else
    txt = {sprintf('Correlation: %3.3f',YpY(x,y)),['Top: ',fname.m{x}],...
      ['Bottom: ',fname.m{y}]};
  end

  % text info for textbox
  if mesh_detected
    txt2 = {sprintf('Correlation: %3.3f',YpY(x,y)),[],'right (1st row) and left (2nd row) hemisphere',['Left: ',...
      spm_file(fname.m{x},'short25') '     Right: ',spm_file(fname.m{y},'short25')],...
      [],'Difference to Sample Mean (red: - green: +)'};
  else
    txt2 = {sprintf('Correlation: %3.3f',YpY(x,y)),[],['Top: ',...
      spm_file(fname.m{x},'short25')],['Bottom: ',spm_file(fname.m{y},'short25')],...
      [],['Displayed slice: ',num2str(round(get(H.mm,'Value'))),' mm'],...
      'Difference to Sample Mean (red: - green: +)'};
  end      

  set(H.text,'String',txt2,'FontSize',FS(6));
  axes('Position',pos.slice);
  x = pos.x;
end
  
if mesh_detected 
  % use indexed 2D-sheet to display surface data as image
  % check surface size to use indexed 2D map
  if (length(data_array(:,x)) == 163842)
    ind = spm_load(fullfile(spm('dir'),'toolbox','cat12','templates_surfaces','fsavg.index2D_256x128.txt'));
    if isscatter
      img = reshape(data_array(ind,x),[256,128]);
    else
      img = [reshape(data_array(ind,x),[256,128]) reshape(data_array(ind,y),[256,128])];
    end
    img = circshift(img,128);
    % alpha overlay
    if isscatter
      img_alpha = reshape(data_array_diff(ind,x),[256,128]);
    else
      img_alpha = [reshape(data_array_diff(ind,x),[256,128]) reshape(data_array_diff(ind,y),[256,128])];
    end
    img_alpha = circshift(img_alpha,128);
  elseif (length(data_array(:,x)) == 327684)
    ind = spm_load(fullfile(spm('dir'),'toolbox','cat12','templates_surfaces','fsavg.index2D_256x128.txt'));
    data_array_x_lh = data_array(1:163842,x);
    data_array_x_rh = data_array(163843:end,x);
    if isscatter
      img_lh = reshape(data_array_x_lh(ind),[256,128]);
      img_rh = reshape(data_array_x_rh(ind),[256,128]);
    else
      data_array_y_lh = data_array(1:163842,y);
      data_array_y_rh = data_array(163843:end,y);
      img_lh = [reshape(data_array_x_lh(ind),[256,128]) reshape(data_array_y_lh(ind),[256,128])];
      img_rh = [reshape(data_array_x_rh(ind),[256,128]) reshape(data_array_y_rh(ind),[256,128])];
    end
    img = [circshift(img_lh,128); img_rh];
    % alpha overlay
    data_array_x_lh = data_array_diff(1:163842,x);
    data_array_x_rh = data_array_diff(163843:end,x);
    if isscatter
      img_lh = reshape(data_array_x_lh(ind),[256,128]);
      img_rh = reshape(data_array_x_rh(ind),[256,128]);
    else
      data_array_y_lh = data_array_diff(1:163842,y);
      data_array_y_rh = data_array_diff(163843:end,y);
      img_lh = [reshape(data_array_x_lh(ind),[256,128]) reshape(data_array_y_lh(ind),[256,128])];
      img_rh = [reshape(data_array_x_rh(ind),[256,128]) reshape(data_array_y_rh(ind),[256,128])];
    end
    img_alpha = [circshift(img_lh,128); img_rh];
  elseif (length(data_array(:,x)) == 32492)
    ind = spm_load(fullfile(spm('dir'),'toolbox','cat12','templates_surfaces_32k','fsavg.index2D_256x128.txt'));
    if isscatter
      img = reshape(data_array(ind,x),[256,128]);
    else
      img = [reshape(data_array(ind,x),[256,128]) reshape(data_array(ind,y),[256,128])];
    end
    img = circshift(img,128);
    % alpha overlay
    if isscatter
      img_alpha = reshape(data_array_diff(ind,x),[256,128]);
    else
      img_alpha = [reshape(data_array_diff(ind,x),[256,128]) reshape(data_array_diff(ind,y),[256,128])];
    end
    img_alpha = circshift(img_alpha,128);
  elseif (length(data_array(:,x)) == 64984)
    ind = spm_load(fullfile(spm('dir'),'toolbox','cat12','templates_surfaces_32k','fsavg.index2D_256x128.txt'));
    data_array_x_lh = data_array(1:32492,x);
    data_array_x_rh = data_array(32493:end,x);
    if isscatter
      img_lh = reshape(data_array_x_lh(ind),[256,128]);
      img_rh = reshape(data_array_x_rh(ind),[256,128]);
    else
      data_array_y_lh = data_array(1:32492,y);
      data_array_y_rh = data_array(32493:end,y);
      img_lh = [reshape(data_array_x_lh(ind),[256,128]) reshape(data_array_y_lh(ind),[256,128])];
      img_rh = [reshape(data_array_x_rh(ind),[256,128]) reshape(data_array_y_rh(ind),[256,128])];
    end
    img = [circshift(img_lh,96); circshift(img_rh,96)];
    % alpha overlay
    data_array_x_lh = data_array_diff(1:32492,x);
    data_array_x_rh = data_array_diff(32493:end,x);
    if isscatter
      img_lh = reshape(data_array_x_lh(ind),[256,128]);
      img_rh = reshape(data_array_x_rh(ind),[256,128]);
    else
      data_array_y_lh = data_array_diff(1:32492,y);
      data_array_y_rh = data_array_diff(32493:end,y);
      img_lh = [reshape(data_array_x_lh(ind),[256,128]) reshape(data_array_y_lh(ind),[256,128])];
      img_rh = [reshape(data_array_x_rh(ind),[256,128]) reshape(data_array_y_rh(ind),[256,128])];
    end
    img_alpha = [circshift(img_lh,96); circshift(img_rh,96)];
  else
    if isscatter
      img = data_array(:,x)';
      % alpha overlay
      img_alpha = data_array_diff(:,x)';
    else
      img = [data_array(:,y) data_array(:,x)]';
      % alpha overlay
      img_alpha = [data_array_diff(:,y) data_array_diff(:,x)]';
    end
  end

  % scale img to 0..64
  mn = min(data_array(:));
  mx = max(data_array(:));
  img = 64*((img - mn)/(mx-mn));
else
  % add slider for colume data
  set(H.mm,'Visible','on');
  set(H.mm_txt,'Visible','on');
  if isscatter
    img = data_array(:,:,x)';
    % alpha overlay
    img_alpha = data_array_diff(:,:,x)';
  else
    img = [data_array(:,:,y) data_array(:,:,x)]';
    % alpha overlay
    img_alpha = [data_array_diff(:,:,y) data_array_diff(:,:,x)]';
  end
end

% correct orientation
img = rot90(img,2);
img_alpha = rot90(img_alpha,2);

% display image with 2nd colorbar (gray)
image(65 + img);
if ~mesh_detected, axis image; end
set(gca,'XTickLabel','','YTickLabel','','TickLength',[0 0]);

% prepare alpha overlays for red and green colors
if alphaval > 0
  % get 2%/98% ranges of difference image
  range = cat_vol_iscaling(img_alpha(:),[0.02 0.98]);

  hold on
  alpha_g = cat(3, zeros(size(img_alpha)), alphaval*ones(size(img_alpha)), zeros(size(img_alpha)));
  alpha_r   = cat(3, alphaval*ones(size(img_alpha)), zeros(size(img_alpha)), zeros(size(img_alpha)));
  hg = image(alpha_g); set(hg, 'AlphaData', img_alpha.*(img_alpha>range(2)),'AlphaDataMapping','scaled')
  if ~mesh_detected, axis image; end
  hr = image(alpha_r); set(hr, 'AlphaData',-img_alpha.*(img_alpha<range(1)),'AlphaDataMapping','scaled')
  if ~mesh_detected, axis image; end
  hold off
end

if mesh_detected
  xlabel('2D surface maps');
end

return
