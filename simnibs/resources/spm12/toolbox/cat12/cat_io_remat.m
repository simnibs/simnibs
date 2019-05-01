function cat_io_remat(P1,Pn)
% ______________________________________________________________________
% Set orientation of Pn by orientation of P1. If more than one P1 input 
% image is used, Pn must have the same number of files.
%
%  cat_io_remat(P1,Pn)
%
%  P1 = image(s) with correct resolution (1 iamge or n images)
%  Pn = goal images for P1.mat (n images) 
% 
% ______________________________________________________________________
% $Revision: 764 $  $Date: 2015-11-17 14:11:53 +0100 (Di, 17 Nov 2015) $

  if ~exist('P1','var') || isempty(P1)
    P1 = spm_select([1 inf],'image','Select correct orientated image');
  else
    P1 = char(P1);
  end
  if numel(P1)==0, return; end
  if ~exist('Pn','var') || isempty(Pn)
    if size(P1,1)==1
      Pn = spm_select(inf,'image','Select uncorrect orientated images');
    else
      Pn = spm_select([size(P1,1) size(P1,1)],'image','Select uncorrect orientated images');
    end  
  else
    Pn = char(Pn);
  end
  if numel(Pn)==0, return; end
  if size(P1,1)>1 && size(P1,1)~=size(P1,1), 
    error('cat_io_repmat:numberP1Pn',...
      sprintf('Number of images P1 and Pn does not fit (n(P1)=%d,n(Pn)=%d).',...
      size(P1,1),size(Pn,1))); %#ok<SPERR>
  end
  
  V1 = spm_vol(P1);
  Vn = spm_vol(Pn);

  for i=1:numel(Vn)
    Y = spm_read_vols(Vn(i));
    if size(P1,1)==1
      Vn(i).mat = V1.mat;
    else
      Vn(i).mat = V1(i).mat;
    end
    spm_write_vol(Vn(i),Y);
  end
end