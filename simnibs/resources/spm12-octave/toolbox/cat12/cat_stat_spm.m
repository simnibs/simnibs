function cat_stat_spm(SPM)
% Workaround to use fsaverage surface as SurfaceID (for displaying results)
% spm_spm is used to estimate the model and the mesh of the 1st file in the model 
% is replaced by the fsaverage brain because this mesh is used for overlaying
% results.
%__________________________________________________________________________
% Christian Gaser
% $Id: cat_stat_spm.m 1239 2017-12-07 14:04:40Z gaser $

if nargin == 0
  P = spm_select([1 Inf],'^SPM\.mat$','Select SPM.mat file(s)');
  for i=1:size(P,1)
    swd = spm_file(P(i,:),'fpath');
    load(fullfile(swd,'SPM.mat'));
    SPM.swd  = swd; 
    cat_stat_spm(SPM);
  end
  
else
  % check for 32k meshes
  if SPM.xY.VY(1).dim(1) == 32492 || SPM.xY.VY(1).dim(1) == 64984
    fsavgDir = fullfile(spm('dir'),'toolbox','cat12','templates_surfaces_32k');
  else
    fsavgDir = fullfile(spm('dir'),'toolbox','cat12','templates_surfaces');
  end
  
  % check that folder exist and number of vertices fits
  if exist(fsavgDir,'dir') == 7 && (SPM.xY.VY(1).dim(1) == 163842 || SPM.xY.VY(1).dim(1) == 327684 || ...
      SPM.xY.VY(1).dim(1) == 655368) || SPM.xY.VY(1).dim(1) == 32492 || SPM.xY.VY(1).dim(1) == 64984
    [pp,ff]   = spm_fileparts(SPM.xY.VY(1).fname);
  
    % find mesh string      
    hemi_ind = [];
    hemi_ind = [hemi_ind strfind(ff,'mesh.')];
    if ~isempty(hemi_ind)
      
      SPM.xY.VY(1).private.private.metadata = struct('name','SurfaceID','value',fullfile(fsavgDir, 'mesh.central.freesurfer.gii'));
      M0 = gifti({fullfile(fsavgDir, 'lh.central.freesurfer.gii'), fullfile(fsavgDir, 'rh.central.freesurfer.gii')});
      G.faces = [M0(1).faces; M0(2).faces+size(M0(1).vertices,1)];
      G.vertices = [M0(1).vertices; M0(2).vertices];
  
      % cerebellar lobes?
      if SPM.xY.VY(1).dim(1) == 655368
        M0 = gifti({fullfile(fsavgDir, 'lc.central.freesurfer.gii'), fullfile(fsavgDir, 'rc.central.freesurfer.gii')});
        G.faces = [G.faces; M0(1).faces+2*size(M0(1).vertices,1); M0(2).faces+3*size(M0(1).vertices,1)];
        G.vertices = [G.vertices; M0(1).vertices; M0(2).vertices];
      end
      
      SPM.xVol.G = G;
      
      % remove memory demanding faces and vertices which are not necessary
      for i=1:length(SPM.xY.VY)
        SPM.xY.VY(i).private.faces = [];
        SPM.xY.VY(i).private.vertices = [];
      end
      
      save(fullfile(SPM.swd,'SPM.mat'),'SPM', '-v7.3');
    else
  
      % find lh|rh string
      hemi_ind = [];
      hemi_ind = [hemi_ind strfind(ff,'lh.')];
      hemi_ind = [hemi_ind strfind(ff,'rh.')];
      hemi = ff(hemi_ind:hemi_ind+1);
      if ~isempty(hemi)
        SPM.xY.VY(1).private.private.metadata = struct('name','SurfaceID','value',fullfile(fsavgDir,[hemi '.central.freesurfer.gii']));
        G = fullfile(fsavgDir,[hemi '.central.freesurfer.gii']);
        SPM.xVol.G = gifti(G);
        
        % remove memory demanding faces and vertices which are not necessary
        for i=1:length(SPM.xY.VY)
          SPM.xY.VY(i).private.faces = [];
          SPM.xY.VY(i).private.vertices = [];
        end
        
        save(fullfile(SPM.swd,'SPM.mat'),'SPM', '-v7.3');
      end
    end
  end
  
  spm_spm(SPM);
end