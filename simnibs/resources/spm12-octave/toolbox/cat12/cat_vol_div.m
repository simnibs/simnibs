function Ydiv = cat_vol_div(Ym,vx_vol)
% ----------------------------------------------------------------------
% Divergence helps to identify all gyri that should not be in the GM, but 
% helps to improve the WM. Divergence estimation is very memory intensive 
% so it is better to limit the resolution.
% ----------------------------------------------------------------------
  if ~exist('vx_vol','var'), vx_vol = repmat(1.5,1,3); end % no reduction
  Ym = single(Ym); 
  [D,I] = cat_vbdist(single(~isnan(Ym))); Ym = Ym(I); % replace nan
  clear D I
  [Ymr,resT2] = cat_vol_resize(Ym,'reduceV',vx_vol,min(1.5,vx_vol*3),32); % don't forget small animals...
  clear Ym
  [gx,gy,gz]  = cat_vol_gradient3(max(1/3,Ymr)); clear Ymr
  
  % divergence function was too memory demanding for some systems
  [px,junk,junk] = cat_vol_gradient3(gx./vx_vol(1)); clear gx junk
  [junk,qy,junk] = cat_vol_gradient3(gy./vx_vol(2)); clear gy junk
  [junk,junk,rz] = cat_vol_gradient3(gz./vx_vol(3)); clear gz junk
  Ydivr = single(px) + single(qy) + single(rz); clear px qy rz
  Ydiv  = cat_vol_resize(smooth3(Ydivr),'dereduceV',resT2); 
return