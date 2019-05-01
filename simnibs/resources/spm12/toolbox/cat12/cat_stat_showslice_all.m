function cat_stat_showslice_all(vargin)
%cat_stat_showslice_all	show 1 slice of all images
%
% FORMAT cat_stat_showslice_all
%
% slice has to be choosen in mm
%__________________________________________________________________________
% Christian Gaser
% $Id: cat_stat_showslice_all.m 1227 2017-11-29 08:20:32Z gaser $

rev = '$Rev: 1227 $';

if nargin == 1
	P = char(vargin.data_vol);
	scaling = vargin.scale;
	slice_mm = vargin.slice;
	orient = vargin.orient;
end

if nargin < 1
	if strcmp(spm('ver'),'SPM2')
		P = spm_get(Inf,'IMAGE','Select normalized files');
	else
		P = spm_select(Inf,'image','Select images');
	end
	scaling = spm_input('Prop. scaling (e.g. for T1- or modulated images)?',1,'yes|no',[1 0],2);
	slice_mm = spm_input('Slice [mm]?','+1','e',0,1);
	orient = spm_input('Orientation',1,'axial|coronal|sagittal',[3 2 1],1);
end

V = spm_vol(deblank(P));
n = size(P,1);

hold = 1;

% voxelsize and origin
vx =  sqrt(sum(V(1).mat(1:3,1:3).^2));
if det(V(1).mat(1:3,1:3))<0, vx(1) = -vx(1); end
Orig = V(1).mat\[0 0 0 1]';

% range
range = ([1 V(1).dim(orient)] - Orig(orient))*vx(orient);

% calculate slice from mm to voxel
sl = slice_mm/vx(orient)+Orig(orient);
while (sl < 1) || (sl > V(1).dim(orient))
	slice_mm = spm_input(['Slice (in mm) [' num2str(range(1)) '...' num2str(range(2)) ']'],1,'e',0);
	sl = slice_mm/vx(orient)+Orig(orient);
end

M = spm_matrix([0 0 0 0 0 0 1 1 1]);
switch orient
  case 1, dim_array = [V(1).dim(2), V(1).dim(3)]; M(:,1:3) = M(:,[2 3 1]);
  case 2, dim_array = [V(1).dim(3), V(1).dim(1)]; M(:,1:3) = M(:,[3 1 2]);
  case 3, dim_array = [V(1).dim(1), V(1).dim(2)]; M(:,1:3) = M(:,[1 2 3]);
end

M(orient,4) = sl;

% global scaling
if scaling
	gm=zeros(size(V,1),1);
	disp('Calculating globals...');
	for i=1:size(V,1), gm(i) = spm_global(V(i)); end
	gm_all = mean(gm);
	for i=1:n
		V(i).pinfo(1:2,:) = gm_all*V(i).pinfo(1:2,:)/gm(i);
	end
end

Y = zeros([dim_array,n]);

%-Start progress plot
%-----------------------------------------------------------------------
spm_progress_bar('Init',n,'volumes completed');
for i=1:n
	d = spm_slice_vol(V(i),M,dim_array,[hold,NaN]);
    if orient == 2
      d = flipud(rot90(d));
    end
	Y(:,:,i) = d;
	spm_progress_bar('Set',i);
end

spm_progress_bar('Clear')

Fgraph = spm_figure('GetWin','Graphics');
FS      = spm('FontSizes');
figure(Fgraph);
spm_figure('Clear',Fgraph);
WIN    = get(gcf,'Position');
WIN    = WIN(3)/WIN(4);
WIN    = WIN/(dim_array(1)/dim_array(2));
sizex  = round(sqrt(n*WIN));
sizey  = round(n/sizex);

while sizex * sizey < n, sizex = sizex + 1; end
if sizex * (sizey-1) >= n, sizey = sizey - 1; end

img = zeros(sizex*dim_array(1),sizey*dim_array(2));

for i = 1:sizex
   for j = 1:sizey
        k = (sizex-i) + sizex*(j-1);
	if k < n
      img((i-1)*dim_array(1)+1:(i)*dim_array(1),(j-1)*dim_array(2)+1:(j)*dim_array(2)) = fliplr(Y(:,:,(k+1)));
	end
  end
end

axes('Position',[0 0 1 0.96])

img = (rot90(img,3));

imagesc(img);

% print filenames

fs = FS(12);
if n>20, fs = FS(8); end
if n>60, fs = FS(6);  end

[tmp, names] = spm_str_manip(char(V.fname),'C');
if isempty(names), names = struct('s',tmp,'m','','e',''); end % case of one images or identical names

fprintf('Compressed filenames: %s\n',tmp);

for i = 1:sizex
  for j = 1:sizey
    k = (sizex-i) + sizex*(j-1);
    if k < n && k >= 0
      text(round(sizex*dim_array(1)-((i-1)*dim_array(1)+(i)*dim_array(1))/2),(j-1)*dim_array(2)+fs+2,names.m{k+1},...
        'FontSize',fs,'Color','r','HorizontalAlignment','center');
    end
  end
end

axis off;
title(sprintf('Slice %dmm: %s*%s',slice_mm,names.s, names.e),'FontSize',FS(10),'FontWeight','Bold');

return
