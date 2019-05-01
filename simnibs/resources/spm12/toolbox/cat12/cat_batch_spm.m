function cat_batch_spm(batchname)
% wrapper for using spm8 batch mode (see cat_batch_cat12.sh)
%_______________________________________________________________________
% $Id: cat_batch_spm.m 767 2015-11-17 16:00:13Z gaser $

if nargin < 1
	fprintf('Syntax: cat_batch_spm(batchname)\n');
	exit
end

spm_get_defaults;
global defaults

if ~exist(batchname,'file')
	fprintf('Batchfile %s not found\n',batchname);
	exit
end

eval(batchname)

if ~exist('matlabbatch','var')
	fprintf('Batchfile %s did not returned variable matlabbatch.\n', batchname);
	exit
end

warning off
spm_jobman('initcfg');
spm_jobman('run',matlabbatch);

exit
