function out = cat_vol_nonlin_coreg_multi_run(job)
% Call cat_vol_nonlin_coreg for multiple subjects
%
% Christian Gaser
% $Id: cat_vol_nonlin_coreg_multi_run.m 1229 2017-11-29 21:03:22Z gaser $

global vox reg bb

warning off;

% use some options from GUI or default file
vox  = job.vox;
reg  = job.reg;
bb   = job.bb;

inputs = cell(3, 1);

m = numel(job.other);
out.ofiles = cell(m,1);
other = cell(m,1);

for j=1:m
  [pth,nam,ext,num] = spm_fileparts(job.other{j});
  out.ofiles{j} = fullfile(pth,['w', nam, ext, num]);
  other{j} = job.other{j};
end
inputs{1} = job.ref;
inputs{2} = job.source;
inputs{3} = other;

spm_jobman('run',{'cat_vol_nonlin_coreg.m'},inputs{:});

warning on;
