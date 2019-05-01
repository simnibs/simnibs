function cat_install_atlases
% Add Dartel atlas labels to spm atlas folder
%
%_______________________________________________________________________
% Christian Gaser
% $Id: cat_install_atlases.m 1133 2017-05-23 08:55:55Z gaser $

spm_dir = spm('dir');
atlas_dir = fullfile(spm_dir,'atlas');
[ST, RS] = mkdir(atlas_dir);

if ST
  [xml_files, n] = cat_vol_findfiles(fullfile(spm_dir,'toolbox','cat12','templates_1.50mm'), 'labels_dartel_*');
  for i = 1:n
    xml_file = deblank(xml_files{i});
    [pth,nam] = spm_fileparts(xml_file);
    atlas_file = fullfile(pth,[nam(15:end) '.nii']);
    try
      copyfile(atlas_file,atlas_dir);
      copyfile(xml_file,atlas_dir);
      fprintf('Install %s\n',atlas_file);
    catch
      disp('Writing error: Please check file permissions.');
    end
  end
else
  error(RS);
end

% this is maybe not enough, to update the file in SPM functions
% you may need to remove the old files and finish SPM, update and restart SPM 
spm_atlas('list','installed','-refresh');

fprintf('\nUse atlas function in SPM Results or context menu in orthogonal view (via right mouse button): Display|Labels\n');