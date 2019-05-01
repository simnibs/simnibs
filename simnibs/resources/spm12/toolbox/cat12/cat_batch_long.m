function cat_batch_long(namefile,output_surface,cat_defaults)
% wrapper for using batch mode (see cat_batch_long.sh)
%
% namefile      - array of file names
% cat_defaults  - use this default file instead of cat_defaults.m
%_______________________________________________________________________
% $Id: cat_batch_long.m 1270 2018-02-07 13:32:14Z gaser $

if nargin < 1
	fprintf('Syntax: cat_batch_long(namefile)\n');
	exit
end

if nargin < 2
  output_surface = 1;
end

fid = fopen(namefile,'r');
names = textscan(fid,'%s');
names = names{:};
fclose(fid);

n = length(names);

if n == 0, error('No file found in %s.\n',namefile); end

global defaults cat matlabbatch

spm_get_defaults;

if nargin < 3
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


matlabbatch{1}.spm.tools.cat.tools.long.subj.mov = cell(n,1);
for i=1:n
  matlabbatch{1}.spm.tools.cat.tools.long.subj.mov{i} = names{i};
end

matlabbatch{1}.spm.tools.cat.tools.long.modulate = 1;
matlabbatch{1}.spm.tools.cat.tools.long.warps = 0;

if output_surface
  matlabbatch{1}.spm.tools.cat.tools.long.output.surface = 1;
end

% always deselect print option
matlabbatch{1}.spm.tools.cat.tools.long.extopts.print = 0;

warning off
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

exit
