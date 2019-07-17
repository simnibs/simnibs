%-----------------------------------------------------------------------
% Job for longitudinal batch
% Christian Gaser
% $Id: cat_long_main.m 1114 2017-03-02 10:46:01Z gaser $
%-----------------------------------------------------------------------

global opts extopts output modulate dartel warps

warning('off','MATLAB:DELETE:FileNotFound');
matlabbatch{1}.spm.tools.cat.tools.series.data = '<UNDEFINED>';

% use some options from gui or default file
for j=2:3
  if exist('opts','var')
    matlabbatch{j}.spm.tools.cat.estwrite.opts = opts;
  end
  if exist('extopts','var')
    matlabbatch{j}.spm.tools.cat.estwrite.extopts = extopts;
  end
  if exist('output','var')
    matlabbatch{j}.spm.tools.cat.estwrite.output = output;
  end
  matlabbatch{j}.spm.tools.cat.estwrite.nproc = 0;
end

% modulation option for applying deformations
if modulate
  matlabbatch{4}.spm.tools.cat.tools.defs.modulate = modulate;
end

% dartel export option
matlabbatch{2}.spm.tools.cat.estwrite.output.GM.dartel = dartel;
matlabbatch{2}.spm.tools.cat.estwrite.output.WM.dartel = dartel;
matlabbatch{3}.spm.tools.cat.estwrite.output.GM.dartel = dartel;
matlabbatch{3}.spm.tools.cat.estwrite.output.WM.dartel = dartel;

matlabbatch{1}.spm.tools.cat.tools.series.bparam = 1000000;
matlabbatch{2}.spm.tools.cat.estwrite.data(1) = cfg_dep('Longitudinal Rigid Registration: Midpoint Average', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','avg', '()',{':'}));
matlabbatch{2}.spm.tools.cat.estwrite.nproc = 0;
matlabbatch{2}.spm.tools.cat.estwrite.output.ROI = 0;
matlabbatch{2}.spm.tools.cat.estwrite.output.surface = 0;
matlabbatch{2}.spm.tools.cat.estwrite.output.GM.mod = 0;
matlabbatch{2}.spm.tools.cat.estwrite.output.WM.mod = 0;
matlabbatch{2}.spm.tools.cat.estwrite.output.bias.warped = 1;
matlabbatch{2}.spm.tools.cat.estwrite.output.warps = [1 0];
matlabbatch{3}.spm.tools.cat.estwrite.data(1) = cfg_dep('Longitudinal Rigid Registration: Realigned images', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','rimg', '()',{':'}));
matlabbatch{3}.spm.tools.cat.estwrite.nproc = 0;
matlabbatch{3}.spm.tools.cat.estwrite.output.GM.native = 1;
matlabbatch{3}.spm.tools.cat.estwrite.output.GM.mod = 0;
matlabbatch{3}.spm.tools.cat.estwrite.output.WM.native = 1;
matlabbatch{3}.spm.tools.cat.estwrite.output.WM.mod = 0;
matlabbatch{3}.spm.tools.cat.estwrite.output.bias.warped = 0;
matlabbatch{4}.spm.tools.cat.tools.defs.field1(1) = cfg_dep('CAT12: Segmentation: Deformation Field', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','fordef', '()',{':'}));
matlabbatch{4}.spm.tools.cat.tools.defs.images(1) = cfg_dep('CAT12: Segmentation: p1 Images', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','tiss', '()',{1}, '.','p', '()',{':'}));
matlabbatch{4}.spm.tools.cat.tools.defs.images(2) = cfg_dep('CAT12: Segmentation: p2 Images', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','tiss', '()',{2}, '.','p', '()',{':'}));
matlabbatch{4}.spm.tools.cat.tools.defs.interp = 1;
matlabbatch{5}.cfg_basicio.file_dir.file_ops.file_move.files(2) = cfg_dep('CAT12: Segmentation: p1 Images', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','tiss', '()',{1}, '.','p', '()',{':'}));
matlabbatch{5}.cfg_basicio.file_dir.file_ops.file_move.files(1) = cfg_dep('CAT12: Segmentation: p2 Images', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','tiss', '()',{2}, '.','p', '()',{':'}));

% save deformations
if ~warps
  matlabbatch{5}.cfg_basicio.file_dir.file_ops.file_move.files(3) = cfg_dep('CAT12: Segmentation: Deformation Field', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','fordef', '()',{':'}));
end

matlabbatch{5}.cfg_basicio.file_dir.file_ops.file_move.action.delete = false;

