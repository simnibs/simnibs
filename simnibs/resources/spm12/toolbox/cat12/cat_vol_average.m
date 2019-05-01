function PO = cat_vol_average(P,fname,PT,dt,nr,mask)
% ______________________________________________________________________
% Creates median images of a set of files P or volumes V with the same
% image properties.
%
%   VO = cat_vol_average(V[,fname,PT,dt,nr])
%   VO = cat_vol_average(P[,fname,PT,dt,nr])
%
%   P     = char or cell array with filenames
%   PT    = template image with another resolution (for interpolation)
%   V     = SPM volume structure
%   VO    = output volume
%   dt    = for [median,mean,std] images with spm nii datatype
%           (0=no image, 2=uint8, 4=int16, 16=single, ...)
%   nr    = use number in fname
%   fname = name of the outputfile, add median/mean/std automaticly
%            mypath/myname.nii > mypath/median003_myname.nii ...
%   mask  = value for masking
% ______________________________________________________________________
% Robert Dahnke 2013_08
% Structural Brain Mapping Group
% University Jena
%  
% $Id: cat_vol_average.m 764 2015-11-17 13:11:53Z gaser $
% ______________________________________________________________________
  
  if isstruct(P)
    V = P; clear P;
    for fi = 1:numel(V), P{fi} = V(fi).fname; end
  end
  P = char(P);
  if isempty(P), VO=struct(); return; end
  
  if ~exist('nr','var') || isempty(nr)
    nr = 0;
  end
  if nr==0
    nstr = '';
  else
    nstr=num2str(size(P,1),'%3.0f');
  end
  
  if ~exist('fname','var') || isempty(fname)
    [pp,ff,ee] = spm_fileparts(P(1,:));
    fname2{1} = fullfile(pp,['median' nstr '_'  ff ee]);
    fname2{2} = fullfile(pp,['mean'   nstr '_'  ff ee]);
    fname2{3} = fullfile(pp,['std'    nstr '_'  ff ee]);
  else
    if iscell(fname)
      fname2 = fname; 
    else
      [pp,ff,ee] = spm_fileparts(fname);
      fname2{1} = fullfile(pp,['median' nstr '_' ff ee]);
      fname2{2} = fullfile(pp,['mean'   nstr '_' ff ee]);
      fname2{3} = fullfile(pp,['std'    nstr '_' ff ee]);      
    end
  end

  if ~exist('dt','var') || isempty(dt)
    % median, mean, sd
    dt = [16 16 16];
  end
  if size(P,1)<2,
    warning('MATLAB:cat_vol_average:input','WARNING:cat_vol_average:to small input (n=1)!\n');
    VO=struct(); return; 
  end
  if size(P,1)<3,
    warning('MATLAB:cat_vol_average:input','WARNING:cat_vol_average:to small input for median and std (n=2)!\n');
    dt(1)=0; dt(1)=0;
  end
  
  
  if numel(dt)==1; dt=repmat(dt,1,3); end
  
  if exist('PT','var') && ~isempty(PT) %&& exist(PT,'file')
  % reslicing
    if iscell(PT) && size(PT,1)<size(PT,2), PT=PT'; end
    PT = char(PT);
    if ~exist(PT,'file')
      error('ERROR:cat_vol_median:PT','ERROR:cat_vol_median:PT-file ''%s'' does not exist.',PT);
    end
    
    para.reslice.mask   = 0;
    para.reslice.mean   = 0;
    para.reslice.interp = 1;
    para.reslice.which  = 1;
    para.reslice.wrap   = [0 0 0];
    para.reslice.prefix = 'median_r';
    
    spm_reslice([cellstr(PT);cellstr(P)],para.reslice);
    
    PR=cell(size(P,1),1);
    for fi = 1:size(P,1)
      [pp,ff,ee] = spm_fileparts(P(fi,:));
      PR{fi} = fullfile(pp,[para.reslice.prefix  ff ee]);
    end
    P = char(PR); clear PR;
    cat_stat_spm_reprint; 
  end
  
  if ~exist('mask','var'), mask=0; end
  
  PO = cell(1,3);
  for i=1:3
    if dt(i)>0
      % median

      
      V  = spm_vol(P); if exist(fname2{i},'file'), delete(fname2{i}); end
      VO1 = V(1); VO1.fname = fname2{i}; VO1.dt(1) = 16;
      switch i
        case 1, VO1.descript = sprintf('median image of %s scans',size(P,1));
        case 2, VO2.descript = sprintf('mean image of %s scans',size(P,1));
        case 3, VO3.descript = sprintf('std image of %s scans',size(P,1));
      end
      
      VO1 = spm_create_vol(VO1);
      Y  = zeros([VO1.dim(1:2),1,size(P,1)],'single');
      for p=1:V(1).dim(3)
        for fi = 1:size(P,1), 
           Y(:,:,1,fi) = single(spm_slice_vol(V(fi),spm_matrix([0 0 p]),V(fi).dim(1:2),0));
        end
        switch i
          case 1, YO = cat_stat_nanmedian(Y,4);
          case 2, YO = cat_stat_nanmean(Y,4);
          case 3, YO = cat_stat_nanstd(Y,4);
        end
        if mask
          YM = cat_stat_nansum(Y>0,4)>=(mask*size(P,1));
          YO = YO .* YM;
        end
        VO1 = spm_write_plane(VO1,YO,p);
      end
      
      if exist('dt','var')
        Y = spm_read_vols(VO1);
        VO1.dt(1) = dt(i); 
        if dt(i)==2 || dt(i)==4
          VO1.pinfo(1) = round(max(Y(:))/(16^dt(i)-1));
        end
        VO1 = rmfield(VO1,'private');
        if exist(fname2{i},'file'), delete(fname2{i}); end
        
        spm_write_vol(VO1,Y);
      end  
      PO{i}=fname2{i};
    end
  end
 
  if exist('PT','file')
    for fi = 1:numel(P)
      delete(P(fi,:))
    end
  end
end
function cat_stat_spm_reprint(str,lines)
  if ~exist('str','var'), str = ''; end
  if ~exist('lines','var'), lines=3; end
  if lines>0
    fprintf(sprintf('%s',repmat('\b',1,lines*73+1)));
  else
    fprintf(sprintf('%s',repmat('\b',1,-lines)));
  end
  fprintf(str);
end
