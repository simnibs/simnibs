function cat_batch_cat(namefile,cat_defaults)
% wrapper for using batch mode (see cat_batch_cat.sh)
%
% namefile      - array of file names
% cat_defaults  - use this default file instead of cat_defaults.m
%
%_______________________________________________________________________
% $Id: cat_batch_cat.m 1270 2018-02-07 13:32:14Z gaser $

 %#ok<*TRYNC>
 
if nargin < 1
	fprintf('Syntax: cat_batch_cat(namefile)\n');
	return
end

[t,pid] = system('echo $$');
fprintf('cat_batch_cat: \n  PID = %s\n\n',pid);

global defaults cat matlabbatch %#ok<NUSED>

spm_get_defaults;

if nargin < 2
    cat_get_defaults;
else
    if isempty(cat_defaults)
        cat_get_defaults;
    else
        fprintf('Use defaults in %s.\n',cat_defaults);
        [pp, name] = spm_fileparts(cat_defaults);
        clear cat_defaults
        oldpath = pwd;
        cd(pp)
        eval(name);
        cd(oldpath)
    end
end

fid = fopen(namefile,'r');
names = textscan(fid,'%s');
names = names{:};
fclose(fid);
n = length(names);

if n == 0, error(sprintf('No file found in %s.\n',namefile)); end %#ok<SPERR>

matlabbatch{1}.spm.tools.cat.estwrite = cat;
matlabbatch{1}.spm.tools.cat.estwrite.data = cellstr(names);

tmp_fields = char('NCstr','BVCstr','regstr','WMHC','WMHCstr','mrf','INV','restype','resval','species','darteltpm','shootingtpm',...
            'cat12atlas','brainmask','T1','pbtres','close_parahipp','scale_cortex','add_parahipp','colormap','verb','ignoreErrors',...
            'expertgui','subfolders','experimental','atlas','LAB','print');
for i=1:size(tmp_fields,1)
  try
    matlabbatch{1}.spm.tools.cat.estwrite.extopts = rmfield(matlabbatch{1}.spm.tools.cat.estwrite.extopts,deblank(tmp_fields(i,:)));
  end
end

tmp_fields = char('atlas','te','pc','WMH','ROI','TPMC','label','CSF','WM','GM','las','bias');
for i=1:size(tmp_fields,1)
  try
    matlabbatch{1}.spm.tools.cat.estwrite.output = rmfield(matlabbatch{1}.spm.tools.cat.estwrite.output,deblank(tmp_fields(i,:)));
  end
end

tmp_fields = char('ngaus','warpreg','biasreg','biasfwhm','samp');
for i=1:size(tmp_fields,1)
  try
    matlabbatch{1}.spm.tools.cat.estwrite.opts = rmfield(matlabbatch{1}.spm.tools.cat.estwrite.opts,deblank(tmp_fields(i,:)));
  end
end

tmp_fields  = char('mod','native','warped','dartel');
tmp_map = char('GM','WM','CSF','bias','las','WMH','label','jacobian');
for i=1:size(tmp_map,1)
  for j=1:size(tmp_fields,1)  
    if isfield(matlabbatch{1}.spm.tools.cat.estwrite.output,(deblank(tmp_map(i,:))))
      if isfield(matlabbatch{1}.spm.tools.cat.estwrite.output.(deblank(tmp_map(i,:))),deblank(tmp_fields(j,:)))
        matlabbatch{1}.spm.tools.cat.estwrite.output.(deblank(tmp_map(i,:))) = rmfield(matlabbatch{1}.spm.tools.cat.estwrite.output.(deblank(tmp_map(i,:))),deblank(tmp_fields(j,:)));
      end
    end
  end
end

% deselect multi-threading for batch
matlabbatch{1}.spm.tools.cat.estwrite.nproc = 0;

try
  spm_jobman('initcfg');
  spm_jobman('run',matlabbatch);
catch %#ok<CTCH> % catch with lasterror is necessary for old matlab versions
  caterr = lasterror;  %#ok<LERR>
  sprintf('\n%s\nCAT Preprocessing error: %s:\n%s\n', repmat('-',1,72),caterr.identifier,caterr.message,repmat('-',1,72));
  for si=1:numel(caterr.stack), cat_io_cprintf('err',sprintf('%5d - %s\n',caterr.stack(si).line,caterr.stack(si).name)); end;
  cat_io_cprintf('err',sprintf('%s\\n',repmat('-',1,72)));  
  error('Batch failed.');
end

spm_unlink(char(namefile))

warning off
exit
