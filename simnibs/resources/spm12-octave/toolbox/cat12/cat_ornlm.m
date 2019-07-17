function out = cat_ornlm(in, v, f, h)
% FORMAT out = cat_ornlm(in, v, f, h)
% 
% Optimized Blockwise Non Local Means Denoising Filter
%
% v - size of search volume (M in paper)
% f - size of neighborhood (d in paper)
% h - smoothing parameter
%
%                          Details on ONLM filter                        
% ***************************************************************************
%  The ONLM filter is described in:                                       
%                                                                         
%  P. Coupe, P. Yger, S. Prima, P. Hellier, C. Kervrann, C. Barillot.     
%  An Optimized Blockwise Non Local Means Denoising Filter for 3D Magnetic
%  Resonance Images. IEEE Transactions on Medical Imaging, 27(4):425-441, 
%  April 2008                                                             
% ***************************************************************************
%
% Christian Gaser
% $Id: cat_ornlm.m 772 2015-11-18 11:02:01Z gaser $

rev = '$Rev: 772 $';

disp('Compiling cat_ornlm.c')

pth = fileparts(which(mfilename));
p_path = pwd;
cd(pth);
mex -O cat_ornlm.c ornlm.c 
cd(p_path);

out = cat_ornlm(in, v, f, h);

return
