function cat_vol_slice_overlay_rgb_ui
%__________________________________________________________________________
% Christian Gaser
% $Id: cat_vol_slice_overlay_rgb_ui.m 1209 2017-11-07 16:20:27Z gaser $

clear global SO
global SO

% ---------------------------------------------------------------------------------------
% image array of max. 3 images
% ---------------------------------------------------------------------------------------
name = char(fullfile(spm('dir'),'toolbox','cat12','templates_1.50mm/cobra.nii'),...
            fullfile(spm('dir'),'toolbox','cat12','templates_1.50mm/mori.nii'),...
            fullfile(spm('dir'),'toolbox','cat12','templates_1.50mm/cat.nii'));

range = [0 1]; % this should be adapted to the image range
logP = 0;      % option to use log-scaled colorbars if the input is a log-transformed p-map

% ---------------------------------------------------------------------------------------
% underlying image
% ---------------------------------------------------------------------------------------
SO.img(1).vol = spm_vol(fullfile(spm('dir'),'toolbox','cat12','templates_1.50mm/Template_T1_IXI555_MNI152_GS.nii'));
SO.img(1).prop = 1;
SO.img(1).cmap = gray;
SO.img(1).range = [0.2 1];  % image range have to be adapted

% ---------------------------------------------------------------------------------------
% selection of slices
% ---------------------------------------------------------------------------------------
% this should be adapted
slices    = [{-42:3:52}, {-70:3:30}];
transform = str2mat('axial','coronal');
sl_name   = str2mat('axial: -42:3:59','coronal: -70:3:34');

% Comment this out if you don't wish slice labels
SO.labels=[];

% ---------------------------------------------------------------------------------------
% ---------------------------------------------------------------------------------------
ind = 1;
nm = deblank(name(ind,:));
range = range(ind,:);
logP = logP(ind);

n = size(sl_name,1);
str = deblank(sl_name(1,:));
for i = 2:n, str = [str '|' deblank(sl_name(i,:))]; end
ind = spm_input('Select slices',1,'m',str);
transform = deblank(transform(ind,:));
slices = slices{ind};


for i=(1:size(name,1))+1
    SO.img(i).background = 0;
	SO.img(i).vol = spm_vol(deblank(name(i-1,:)));
	SO.img(i).prop = 1;

	cmap = zeros(64,3);
	cmap(:,i-1) = (1:64)'/64;	
	SO.img(i).cmap = cmap;

	if range(1)==range(2)
		[mn mx] = cg_max(SO.img(i).vol);
		SO.img(i).range = [mn mx];
	else SO.img(i).range = range; end
	SO.img(i).func = 'i1(i1==0)=NaN;';

	if range(1) > 0
		SO.img(i).outofrange = {0,size(SO.img(i).cmap,1)};
	else
		SO.img(i).outofrange = {1,0};
	end
end

SO.transform = transform;
SO.slices = slices;
SO.cbar=[(1:size(name,1))+1];

n_images = length(slices) + length(SO.cbar);
xy = get_xy(n_images);

n = size(xy,1);
xy_name = num2str(xy);
str = deblank(xy_name(1,:));
for i = 2:n, str = [str '|' deblank(xy_name(i,:))]; end
ind = spm_input('Select xy',1,'m',str);
xy = xy(ind,:);


SO.xslices = xy(:,1);
switch lower(SO.transform)
	case 'sagittal'
		dim = xy.*SO.img(1).vol.dim(2:3);
	case 'coronal'
		dim = xy.*SO.img(1).vol.dim([1 3]);
	case 'axial'
		dim = xy.*SO.img(1).vol.dim(1:2);
end
screensize = get(0,'screensize');

scale = screensize(3:4)./dim;
% scale image only if its larger than screensize
if min(scale) < 1
	fig_size = min(scale)*dim*0.975;
else
	fig_size = dim;
end

h = figure(12);
set(h,...
	'Position',[1 1 fig_size],...
	'MenuBar','none',...
	'Resize','off',...
	'PaperType','A4',...
	'PaperUnits','normalized',...
	'PaperPositionMode','auto',...
	'Visible','off');

SO.figure = h;
SO.area.units='pixels';

slice_overlay

% save image
saving = spm_input('Save png images?','+1','yes|no',[1 0],2);
if saving
	[pt,nm] = fileparts(nm);
	imaname = spm_input('Filename','+1','s',[nm '_' lower(transform) '.png']);
	slice_overlay('print',imaname,'print -dpng -painters -noui')
	fprintf('Image %s saved.\n',imaname);
end

return

function xy=get_xy(n)

nn = round(n^0.4);
if n>8, x = nn:round(n/nn); else x = 1:n; end
xy=[];
for i=1:length(x)
	y = round(n/x(i));
	% check whether y is to small
	while y*x(i)<n, y = y + 1; end
	if i>2
		if y*x(i-1)<n, xy = [xy; [x(i) y]]; end
	else xy = [xy; [x(i) y]]; end
end

% change order of x and y
for i=1:size(xy,2)
	yx = [xy(i,2) xy(i,1)];
	xy = [xy; yx];
end

% remove duplicates
xy = unique(xy,'rows');
return

% --------------------------------------------------------------------------
function s=remove_zeros(s)

pos = length(s);
while pos>1
	if strcmp(s(pos),'0')
		s(pos)='';
		pos = pos-1;
	else break
	end
end
