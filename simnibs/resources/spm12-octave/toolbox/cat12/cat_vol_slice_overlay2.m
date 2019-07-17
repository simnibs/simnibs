function cat_vol_slice_overlay(OV)
%__________________________________________________________________________
% Christian Gaser
% $Id: cat_vol_slice_overlay.m 1234 2017-12-04 10:52:18Z gaser $

clear global SO
global SO

if nargin == 0
    
    imgs = spm_select(2, 'image', 'Select additional overlay image', {fullfile(spm('dir'), 'toolbox', 'cat12', 'templates_1.50mm/Template_T1_IXI555_MNI152_GS.nii')});
    OV = pr_basic_ui(imgs, 0);
    
    % set options
    OV.opacity = 1;
    OV.reference_image = deblank(imgs(1, :));
    OV.reference_range = OV.img(1).range;
    OV.name = imgs(2:end, :);
    OV.cmap = OV.img(2).cmap;
    OV.range = OV.img(2).range;
    OV.slices_str = '';
    
end

FS = spm('FontSizes');

% check filename whether log. scaling was used
OV.logP = zeros(size(OV.name, 1));
for i = 1:size(OV.name, 1)
    if contains(OV.name(i, :), 'logP') || contains(OV.name(i, :), 'log_')
        OV.logP(i) = 1;
    end
end

% check fields of OV structure
fieldnames = char('reference_image', 'reference_range', ...
    'opacity', 'cmap', 'name', 'range', 'logP', 'slices_str', 'transform');
for i = 1:size(fieldnames, 1)
    str = deblank(fieldnames(i, :));
    if ~isfield(OV, str)
        error([str, ' not defined']);
    end
end

cmap_bivariate = [1 - (hot); hot]; % colormap if range(1) < 0

if isfield(OV, 'labels')
    SO.labels = OV.labels;
end

if isfield(OV, 'cbar')
    SO.cbar = OV.cbar;
else
    SO.cbar = 2; % colorbar
end

n = size(OV.name, 1);

str = deblank(OV.name(1, :));
for i = 2:n, str = [str, '|', deblank(OV.name(i, :))]; end

if n > 1
    sel = spm_input('Select image', 1, 'm', str);
else
    sel = 1;
end

nm = deblank(OV.name(sel, :));

% if only one argument is given assume that parameters are the same for all files
if size(OV.range, 1) > 1
    range = OV.range(sel, :);
else
    range = OV.range;
end

if size(OV.logP, 1) > 1
    logP = OV.logP(sel);
else
    logP = OV.logP;
end

img = nm;

n_slice = size(OV.slices_str, 1);
if n_slice > 0
    for i = 1:n_slice
        try
            slices{i} = eval(OV.slices_str(i, :));
        catch
            slices{i} = [];
        end
    end
else
    slices{1} = OV.slices;
end

sl_name = [];
for i = 1:size(OV.transform, 1)
    if n_slice > 0
        sl_name = strvcat(sl_name, [OV.transform(i, :), ': ', OV.slices_str(i, :)]);
    else
        sl_name = strvcat(sl_name, OV.transform(i, :));
    end
end

str_select = deblank(sl_name(1, :));
for i = 2:n_slice, str_select = [str_select, '|', deblank(sl_name(i, :))]; end
ind = spm_input('Select slices', '+1', 'm', str_select);
OV.transform = deblank(OV.transform(ind, :));
slices = slices{ind};

SO.img(1).vol = spm_vol(OV.reference_image);
SO.img(1).prop = 1;
SO.img(1).cmap = gray;
SO.img(1).range = OV.reference_range;
SO.img(1).background = 0;

SO.img(2).vol = spm_vol(img);

SO.img(2).prop = OV.opacity; % transparent overlay
SO.img(2).cmap = OV.cmap; % colormap
if ~isfield(OV, 'func')
    SO.img(2).func = 'i1(i1==0)=NaN;';
else
    SO.img(2).func = OV.func;
end

if ~isfield(OV, 'range')
    [mx, mn] = volmaxmin(OV.img(2).vol);
    SO.img(2).range = spm_input('Intensity range for colormap', '+1', 'e', [mn, mx], 2)';
else
    SO.img(2).range = range;
end

if range(1) >= 0
    SO.img(2).outofrange = {0, size(SO.img(2).cmap, 1)};
else
    SO.img(2).outofrange = {1, 1};
    SO.img(2).cmap = cmap_bivariate;
end

SO.transform = OV.transform;
SO.slices = slices;

if isempty(SO.slices)
    [mx, mn, XYZ, vol] = volmaxmin(SO.img(2).vol);
    
    % threshold map and restrict coordinates
    Q = find(vol >= SO.img(2).range(1) & vol <= SO.img(2).range(2));
    XYZ = XYZ(:, Q);
    vol = vol(Q);
    
    M = SO.img(2).vol.mat;
    XYZmm = M(1:3, :) * [XYZ; ones(1, size(XYZ, 2))];
    orientn = strcmpi(SO.transform, {'sagittal', 'coronal', 'axial'});
    
    XYZ_unique = get_xyz_unique(XYZ, XYZmm, vol);
    SO.slices = XYZ_unique{orientn};
end

n_images = length(SO.slices) + length(SO.cbar);
xy = get_xy(n_images);

n = size(xy, 1);
xy_name = num2str(xy);
str = deblank(xy_name(1, :));
for i = 2:n, str = [str, '|', deblank(xy_name(i, :))]; end

if ~isfield(OV, 'xy')
    indxy = spm_input('Select number of columns/rows', '+1', 'm', str);
    xy = xy(indxy, :);
else
    xy = OV.xy;
end

% prepare overview of slices
V = SO.img(1).vol;
ref_vol = spm_read_vols(V);
ref_vol = 64 * (ref_vol - SO.img(1).range(1)) / (SO.img(1).range(2) - SO.img(1).range(1));
vx = sqrt(sum(V.mat(1:3, 1:3).^2));
Orig = round(V.mat \ [0, 0, 0, 1]');

h0 = figure(11);
clf
axes('Position', [0, 0, 1, 1]);

hold on
dim = SO.img(1).vol.dim(1:3);
switch lower(OV.transform)
    case 'sagittal'
        ref_img = ref_vol(:, :, Orig(3))';
        slices_vx = SO.slices / vx(1) + Orig(1);
        image(fliplr(ref_img))
        for i = slices_vx
            h = line([i, i], [1, dim(2)]);
            set(h, 'Color', 'r')
        end
    case 'coronal'
        ref_img = squeeze(ref_vol(Orig(1), :, :))';
        slices_vx = SO.slices / vx(2) + Orig(2);
        image(ref_img)
        for i = slices_vx
            h = line([i, i], [1, dim(3)]);
            set(h, 'Color', 'r')
        end
    case 'axial'
        ref_img = squeeze(ref_vol(Orig(1), :, :))';
        slices_vx = SO.slices / vx(3) + Orig(3);
        image(ref_img)
        for i = slices_vx
            h = line([1, dim(2)], [i, i]);
            set(h, 'Color', 'r')
        end
end

% get position of graphic figure
pos1 = get(spm_figure('FindWin', 'Graphics'), 'Position');

screensize = get(0, 'screensize');
set(h0, 'Position', [pos1(1), pos1(2) + 0.9 * screensize(4), 2 * size(ref_img, 2), 2 * size(ref_img, 1)], ...
    'MenuBar', 'none', ...
    'Resize', 'off', ...
    'PaperType', 'A4', ...
    'PaperUnits', 'normalized', ...
    'NumberTitle', 'off', ...
    'Name', 'Slices', ...
    'PaperPositionMode', 'auto');

hold off
axis off xy image
colormap(gray)

SO.xslices = xy(:, 1);
switch lower(OV.transform)
    case 'sagittal'
        dim = xy .* SO.img(1).vol.dim(2:3);
    case 'coronal'
        dim = xy .* SO.img(1).vol.dim([1, 3]);
    case 'axial'
        dim = xy .* SO.img(1).vol.dim(1:2);
end
screensize = get(0, 'screensize');

% use double size
dim = 2 * dim;

scale = screensize(3:4) ./ dim;
% scale image only if its larger than screensize
if min(scale) < 1
    fig_size = min(scale) * dim * 0.975;
else
    fig_size = dim;
end

[pt, nm] = spm_fileparts(img);

h = figure(21);
set(h, ...
    'Position', [pos1(1), pos1(2), fig_size], ...
    'MenuBar', 'none', ...
    'Resize', 'off', ...
    'PaperType', 'A4', ...
    'PaperUnits', 'normalized', ...
    'PaperPositionMode', 'auto', ...
    'NumberTitle', 'off', ...
    'Name', nm, ...
    'Visible', 'off');

SO.figure = h;
SO.area.units = 'pixels';

slice_overlay

% change labels of colorbar for log-scale
H = gca;
if (SO.cbar == 2) && logP
    YTick = get(H, 'YTick');
    mn = floor(min(YTick));
    mx = ceil(max(YTick));
    
    % only allow integer values
    values = floor(mn:mx);
    pos = get(get(gca, 'YLabel'), 'position');
    pos(1) = 2.5;
    
    set(H, 'YTick', values);
    YTick = get(H, 'YTick');
    
    YTickLabel = [];
    for i = 1:length(YTick)
        YTickLabel = char(YTickLabel, remove_zeros(sprintf('%.g', 10^(-YTick(i)))));
    end
    set(H, 'YTickLabel', YTickLabel)
    
    set(get(gca, 'YLabel'), 'string', 'p-value', 'position', pos, 'FontSize', FS(14))
    
end

set(H, 'FontSize', 0.8*get(H, 'FontSize'))

% save image
if ~isfield(OV, 'save')
    image_ext = spm_input('Save image file?', '+1', 'no|png|jpg|pdf|tif', char('none', 'png', 'jpeg', 'pdf', 'tiff'), 2);
else
    if isempty(OV.save)
        image_ext = 'none';
    else
        [pp, nn, ee] = spm_fileparts(OV.save);
        image_ext = ee(2:end);
    end
end

if ~strcmp(image_ext, 'none')
    
    [pt, nm] = spm_fileparts(img);
    
    % use shorter ext for jpeg
    if ~isfield(OV, 'save')
        if strcmp(image_ext, 'jpeg')
            imaname = spm_input('Filename', '+1', 's', [nm, '_', lower(OV.transform), '.jpg']);
        else
            imaname = spm_input('Filename', '+1', 's', [nm, '_', lower(OV.transform), '.', image_ext]);
        end
    else
        imaname = OV.save;
    end
    
    % and print
    H = findobj(get(SO.figure, 'Children'), 'flat', 'Type', 'axes');
    set(H, 'Units', 'normalized')
    
    saveas(SO.figure, imaname, image_ext);
    fprintf('Image %s saved.\n', imaname);
    if n_slice > 0
        imaname = [lower(OV.transform), '_', replace_strings(OV.slices_str(ind, :)), '.', image_ext];
    else
        imaname = [lower(OV.transform), '.', image_ext];
    end
    
    saveas(h0, imaname, image_ext);
    fprintf('Image %s saved.\n', imaname);
end


% --------------------------------------------------------------------------
function xy = get_xy(n)

nn = round(n^0.4);
if n > 8, x = nn:round(n/nn); else x = 1:n; end
xy = [];
for i = 1:length(x)
    y = round(n/x(i));
    % check whether y is to small
    while y * x(i) < n, y = y + 1; end
    if i > 2
        if y * x(i-1) < n, xy = [xy; [x(i), y]]; end
    else xy = [xy; [x(i), y]]; end
end

% change order of x and y
for i = 1:size(xy, 2)
    yx = [xy(i, 2), xy(i, 1)];
    xy = [xy; yx];
end

% remove duplicates
xy = unique(xy, 'rows');
return

% --------------------------------------------------------------------------
function s = remove_zeros(s)

pos = length(s);
while pos > 1
    if strcmp(s(pos), '0')
        s(pos) = '';
        pos = pos - 1;
    else break
    end
end


% --------------------------------------------------------------------------
function s = replace_strings(s)

s = deblank(s);
% replace spaces with "_" and characters like "<" or ">"
s(findstr(s, ' ')) = '_';
s(findstr(s, ':')) = '_';
s = spm_str_manip(s, 'v');

return

% --------------------------------------------------------------------------
function [mx, mn, XYZ, img] = volmaxmin(vol)

if nargout > 2
    XYZ = [];
end
if nargout > 3
    img = [];
end

mx = -Inf; mn = Inf;
for i = 1:vol.dim(3)
    tmp = spm_slice_vol(vol, spm_matrix([0, 0, i]), vol.dim(1:2), [0, NaN]);
    tmp1 = tmp(isfinite(tmp(:)) & (tmp(:) ~= 0));
    if ~isempty(tmp1)
        if nargout > 2
            [Qc, Qr] = find(isfinite(tmp) & (tmp ~= 0));
            if size(Qc, 1)
                XYZ = [XYZ; [Qc, Qr, i * ones(size(Qc))]];
                if nargout > 3
                    img = [img; tmp1];
                end
            end
        end
        mx = max([mx; tmp1]);
        mn = min([mn; tmp1]);
    end
end

if nargout > 2
    XYZ = XYZ';
end

return

% --------------------------------------------------------------------------
function SO = pr_basic_ui(imgs, dispf)
% GUI to request parameters for slice_overlay routine
% FORMAT SO = pr_basic_ui(imgs, dispf)
%
% GUI requests choices while accepting many defaults
%
% imgs  - string or cell array of image names to display
%         (defaults to GUI select if no arguments passed)
% dispf - optional flag: if set, displays overlay (default = 1)
%
% $Id: cat_vol_slice_overlay.m 1234 2017-12-04 10:52:18Z gaser $

if nargin < 1
    imgs = '';
end
if isempty(imgs)
    imgs = spm_select(Inf, 'image', 'Image(s) to display');
end
if ischar(imgs)
    imgs = cellstr(imgs);
end
if nargin < 2
    dispf = 1;
end

clear global SO
global SO

spm_input('!SetNextPos', 1);

% load images
nimgs = length(imgs);

% process names
nchars = 20;
imgns = spm_str_manip(imgs, ['rck', num2str(nchars)]);

SO.transform = deblank(spm_input('Image orientation', '+1', ['Axial|', ...
    ' Coronal|Sagittal'], strvcat('axial', 'coronal', 'sagittal'), ...
    1));
orientn = find(strcmpi(SO.transform, {'sagittal', 'coronal', 'axial'}));

% identify image types
SO.cbar = [];
XYZ_unique = cell(3, 1);
for i = 1:nimgs
    SO.img(i).vol = spm_vol(imgs{i});
    if i == 1
        SO.img(i).cmap = gray;
        [mx, mn] = volmaxmin(SO.img(i).vol);
        SO.img(i).range = [0.15 * mx, 0.9 * mx];
    else
        [mx, mn] = volmaxmin(SO.img(i).vol);
        
        SO.img(i).func = 'i1(i1==0)=NaN;';
        SO.img(i).prop = Inf;
        SO.cbar = [SO.cbar, i];
        SO.img(i).cmap = return_cmap('Colormap:', 'jet');
        SO.img(i).range = spm_input('Img val range for colormap', '+1', 'e', [mn, mx], 2)';
        
        define_slices = spm_input('Slices', '+1', 'm', 'Estimate slices with local maxima|Define slices', [0, 1], 1);
        if ~define_slices
            [mx, mn, XYZ, img] = volmaxmin(SO.img(i).vol);
            % threshold map and restrict coordinates
            Q = find(img >= SO.img(i).range(1) & img <= SO.img(i).range(2));
            XYZ = XYZ(:, Q);
            img = img(Q);
            
            M = SO.img(i).vol.mat;
            XYZmm = M(1:3, :) * [XYZ; ones(1, size(XYZ, 2))];
            
            XYZ_unique = get_xyz_unique(XYZ, XYZmm, img);
        end
        
    end
end


% slices for display
ts = [0, 0, 0, pi / 2, 0, -pi / 2, -1, 1, 1; ...
    0, 0, 0, pi / 2, 0, 0, 1, -1, 1; ...
    0, 0, 0, 0, 0, 0, 1, 1, 1];

V = SO.img(2).vol;
D = V.dim(1:3);
T = spm_matrix(ts(orientn, :)) * V.mat;
vcorners = [1, 1, 1; D(1), 1, 1; 1, D(2), 1; D(1:2), 1; ...
    1, 1, D(3); D(1), 1, D(3); 1, D(2:3); D(1:3)]';
corners = T * [vcorners; ones(1, 8)];

SO.slices = spm_input('Slices to display (mm)', '+1', 'e', XYZ_unique{orientn});
SO.figure = figure(21);

% and do the display
if dispf, slice_overlay; end

return

% --------------------------------------------------------------------------
function cmap = return_cmap(prompt, defmapn)
cmap = [];
while isempty(cmap)
    cmap = slice_overlay('getcmap', spm_input(prompt, '+1', 's', defmapn));
end
return

% --------------------------------------------------------------------------
function XYZ_unique = get_xyz_unique(XYZ, XYZmm, img)

xyz_array = [];

% cluster map
A = spm_clusters(XYZ);
for j = 1:max(A)
    ind = find(A == j);
    xyz = XYZmm(:, ind);
    xyz_array = [xyz_array, xyz(:, img(ind) == max(img(ind)))];
end

% only keep unique coordinates
XYZ_unique = cell(3, 1);
for j = 1:3
    XYZ_unique{j} = unique(xyz_array(j, :));
end

return

