function [Ym,Yb,WMth,Affine] = cat_long_APP(PF,PG,PB,opt)
% Preprocessing for longitudinal pipeline based on the cat_run_job
% APP pipeline.
%
% [Ym,Yb,WMth,Affine] = cat_long_APP(PF,PG,PB,opt)
% 
% PF    .. original image
% PG    .. t1 template for affine registration
% PB    .. initial template mask
%
% opt         .. parameter (see code)
% opt.verb    .. be verbose (0 - no, 1 - less, 2 - more)
% opt.gcutstr .. strength of skullstripping (0-less,1-strong,def.=0.5)
% opt.vx_vol  .. voxel resolution (def.=ones(1,3)) 
% opt.samp    .. sampling distance of affine registration (def.=3); 
%  
% Ym     .. bias corrected image
% Yb     .. new (smooth) brain mask with a range 0..1 (needs to be thresholded)
% WMth   .. WM threshold of the original image
% Affine .. Affine registration matrix
%
% Call in cat_run_job:
%   [Ym,Yb,WMth] = cat_long_run_APP(job.channel(1).vols{subj},...
%     job.extopts.T1,job.extopts.brainmask)
%
% Display result
%   ds('l2','',vx_vol,Ym,Yb,Yp0,Ym,160)
% ______________________________________________________________________
% Robert Dahnke
% $Id: cat_long_APP.m 1121 2017-03-31 06:54:16Z gaser $

%#ok<*WNOFF,*WNON>

  VF = spm_vol(char(PF));
  VG = spm_vol(char(PG));
  VB = spm_vol(char(PB));

  % parameter
  if ~exist('opt','var'), opt = struct(); end
  def.verb    = 2;
  def.gcutstr = 0.5;
  def.vx_vol  = ones(1,3); 
  def.samp    = 3; 
  opt = cat_io_checkinopt(opt,def);

  if opt.verb
    stime = cat_io_cmd('APP'); 
  end
  
  % Rescale images so that globals are better conditioned
  VF.pinfo(1:2,:) = VF.pinfo(1:2,:)/spm_global(VF);
  VG.pinfo(1:2,:) = VG.pinfo(1:2,:)/spm_global(VG);

  % initial APP
  Ysrc = single(VF.private.dat(:,:,:)); 
  [Ym,Yt,Ybg,WMth] = cat_run_job_APP_init(Ysrc,opt.vx_vol,struct('verb',opt.verb-1));
  
  % write data to VF
  VF.dt         = [spm_type('UINT8') spm_platform('bigend')];
  VF.dat(:,:,:) = cat_vol_ctype(Ym * 200,'uint8'); 
  VF.pinfo      = repmat([1;0],1,size(Ym,3));
  clear Yt; 

  % smoothing
  resa  = opt.samp*2; % definine smoothing by sample size
  VF   = spm_smoothto8bit(VF,resa);
  VG   = spm_smoothto8bit(VG,resa);
  
  % prepare affine parameter 
  aflags     = struct('sep',opt.samp,'regtype','subj','WG',[],'WF',[],'globnorm',1); 
  aflags.sep = max(aflags.sep,max(sqrt(sum(VG(1).mat(1:3,1:3).^2))));
  aflags.sep = max(aflags.sep,max(sqrt(sum(VF(1).mat(1:3,1:3).^2))));

  % Affine registration
  try 
    warning off
    [Affine, affscale]  = spm_affreg(VG, VF, aflags, eye(4)); 
    warning on
    clear VG 
  catch
    affscale = 0; 
  end
  if affscale>3 || affscale<0.5 % failed registration
    Affine = eye(4); 
  end
  
  % apply (first affine) registration on the default brain mask
  VFa = VF; VFa.mat = Affine * VF.mat; 
  if isfield(VFa,'dat'), VFa = rmfield(VFa,'dat'); end
  [pp,ff] = spm_fileparts(PF); Pbt = fullfile(pp,['brainmask_' ff '.nii']);
  [Vmsk,Yb]   = cat_vol_imcalc([VFa,VB],Pbt,'i2',struct('interp',3,'verb',0));
    
  if opt.verb
    fprintf('%4.0fs\n',etime(clock,stime)); 
  end
end




















