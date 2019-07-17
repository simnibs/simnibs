function [c,g,w]=cat_io_seg2cgw(seg)
% ______________________________________________________________________
% Convert segment image into 3 images for CSF, GM and WM.
%
%   [c,g,w]=cat_io_seg2cgw(seg)
%
% seg:      p0 tissue segment image (0=BG,1=CSF,2=GM,3=WM)
% [c,g,w]:  tissue class map with values from 0 to 1
% ______________________________________________________________________
% Robert Dahnke
% Structural Brain Mapping Group
% University Jena
%
% $Id: cat_io_seg2cgw.m 764 2015-11-17 13:11:53Z gaser $

  %if isa(char), fseg=seg; hseg=spm_vol(seg); seg=spm_read_vols(hseg); end
  
  c =  seg    .* (seg>0 & seg<1)   +   (2-seg).*(seg>=1 & seg<2);
  g = (seg-1) .* (seg>1 & seg<2)   +   (3-seg).*(seg>=2 & seg<3);
  w = (seg-2) .* (seg>2 & seg<3)   +   (4-seg).*(seg>=3 & seg<4);
  
  % for some loops in other functions...
  %{
  for ci=1:3
    c1 = (v1-(ci-1)).* (v1>(ci-1) & v1<ci) + ((ci+1)-v1).*(v1>=ci & v1<(ci+1));
  end
  %}
  
end