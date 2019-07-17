function Yg = cat_vol_grad(Ym,vx_vol)
% ----------------------------------------------------------------------
% gradient map for edge description
% ----------------------------------------------------------------------
  if ~exist('vx_vol','var'), vx_vol=ones(1,3); end
  Ym = single(Ym); 
  [D,I] = cat_vbdist(single(~isnan(Ym))); Ym = Ym(I); % replace nan
  [gx,gy,gz] = cat_vol_gradient3(single(Ym)); 
  Yg = abs(gx./vx_vol(1))+abs(gy./vx_vol(2))+abs(gz./vx_vol(3)); 
  %Yg = Yg ./ (Ym+eps);
return
