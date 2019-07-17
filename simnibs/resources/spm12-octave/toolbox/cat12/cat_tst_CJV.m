function CJV = cat_tst_CJV(P,Pp0)
% ______________________________________________________________________
% Function to estimate the CJV (coefficient of joint variation) in
% images.
% 
%  CJV = cat_tst_CJV(P,Pp0)
%  
%  P    ... set of n images (cellstr or char)
%  Pp0  ... set of 1 (ground truth) or n images (cellstr of char)
%  CJV  ... matrix of n CJV values of each image
%
% ______________________________________________________________________
% Robert Dahnke 
% Structural Brain Mapping Group
% University Jena
%
% $Id: cat_tst_CJV.m 1171 2017-08-31 08:41:59Z gaser $
% ______________________________________________________________________

  if iscell(P)   && size(P,1)  <size(P,2), end% P=char(P);     end
  if iscell(Pp0) && size(Pp0,1)<size(P,2), end%Pp0=char(Pp0); end

  
  % get headers
  V   = spm_vol(char(P));
  Vp0 = spm_vol(char(Pp0));
  
  % get group gt-image
  if numel(Vp0)==1
    Yp0 = spm_read_vols(Vp0);
  end
  
  % get images
  CJV = zeros(size(V));
  if ~isempty(P)
    for i=1:numel(V)
      try
        Y = spm_read_vols(V(i));
        if numel(Vp0)>1
          Yp0 = spm_read_vols(Vp0(i));
        else
          Yp0 = spm_read_vols(Vp0);
        end
        ncls = max(round(Yp0(:))); 
        if ncls==254
          Yp0 = Yp0/ncls*3;
        end

        CJV(i) = ( cat_stat_nanstd(Y(Yp0>2.5))/2 + ...
                   cat_stat_nanstd(Y(Yp0>1.5 & Yp0<2.5))/2 )  / ...
                 ( cat_stat_nanmean(Y(Yp0>2.5)) - ...
                   cat_stat_nanmean(Y(Yp0>1.5 & Yp0<2.5)) );
      catch
        cat_io_cprintf('err',sprintf('Error: %s\n',P(i,:)));
        CJV(i) = nan;
      end
    end
  else
    cat_io_cprintf('err',sprintf('Error: no images\n'));
  end
end