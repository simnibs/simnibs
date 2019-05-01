function cat_io_volctype(varargin)
% ______________________________________________________________________
% Convert datatype of images, to have more space on your harddisk. 
% In example most tissue classifcations can saved as uint8 or uint16 
% rather than double or single. If the image contain negative values
% int8/16 rather than uint8/16 is used. 
% ______________________________________________________________________
%
% Robert Dahnke 
% Structural Brain Mapping Group
% University Jena
% ______________________________________________________________________
% $Id: cat_io_volctype.m 1170 2017-08-30 16:59:05Z dahnke $

%% choose images
  
  if nargin == 0 
      job.data = cellstr(spm_select([1 Inf],'image','Select images')); 
  else
      job = varargin{1};
  end
  def.prefix    = ''; 
  def.force     = 1; 
  def.overwrite = 1;
  job = cat_io_checkinopt(job,def); 

  if ~isfield(job,'data') || isempty(job.data)
     job.data = cellstr(spm_select([1 Inf],'image','Select images')); 
  else
     job.data = cellstr(job.data);
  end
  if isempty(job.data) || isempty(job.data{1}), return; end
  
% choose output format
  spm_clf('Interactive'); 
  if isfield(job,'ctype')
    if ischar(job.ctype)
      ctype = spm_type(job.ctype);
    else
      ctype = job.ctype;
    end  
  else
    V = spm_vol(strrep(job.data{1},',1',''));
    switch V.dt(1)
      case {2,512}, dtype = 1; % uint8
      case {4,256}, dtype = 3; % int8
      otherwise,    dtype = 2; % uint16
    end
    ctype = spm_input('Datatype',1,'(u)int8|(u)int16|single',[2,512,16],dtype);
  end
  if any(ctype==[2 4 256 512])
    if isfield(job,'range')
      range = job.range; 
    else
      range = spm_input('Range','+1','100%|99.99%|%|ss',[100,99.99,-1,0],2);
    end
    
    if range==0
      V = spm_vol(strrep(job.data{1},',1',''));
      Y = spm_read_vols(V);
      cvals = 1/round(single(intmax(spm_type(ctype))) * diff([min(Y(:)),max(Y(:))]));
    elseif range<0 || range>100
      range = min(100,max(eps,spm_input('Range in %:','+1','r',99.99,1)));
      cvals = 0;
    else
      cvals = 0;
    end
  end
  
  % stepsize
  if ~isfield(job,'cvals') && range==0
    job.cvals = spm_input(sprintf('Stepsize (0=auto;min:%4.2f;max:%4.2f):',min(Y(:)),max(Y(:))),'+1','r',0,1);
  else
    job.cvals = cvals; 
  end
  
% choose prefix
  job.prefix = spm_input('Filename prefix (empty=overwrite!)','+1','s',job.prefix,1);
  %job.overwrite = spm_input('Overwrite?','+1','y/n','',job.overwrite);
  
 
  
  %% convert
  for si=1:numel(job.data)
    V = spm_vol(strrep(job.data{si},',1',''));
    Y = spm_read_vols(V); 
 
    
    [pp,ff,ee,dd] = spm_fileparts(V(1).fname);
      

    if range==0
      clim = [min(Y(:)) max(Y(:))]; 
    elseif range>0 && range<100
      clim = cat_vol_iscaling(Y(:),[50-range/2 50+range/2]);         
    end
    
    if round(cvals)~=cvals && ccvals~=0
      switch ctype
        case {2,256}, V(1).pinfo(1) = (clim(2) - clim(1)) / 256;
        case {4,512}, V(1).pinfo(1) = (clim(2) - clim(1)) / 256^2;
      end
    else
       V(1).pinfo(1) = cvals; 
    end
    
    if any(clim<0)
      switch V(1).dt(1) 
        case 2,   ctype = 4;  
        case 512, ctype = 256; 
      end
    else
      switch V(1).dt(1) 
        case 4,   ctype = 2;  
        case 256, ctype = 512; 
      end
    end
    V(1).dt(1) = ctype;

    if isempty(job.prefix) && (job.overwrite==0 || job.overwrite=='n')
      prefix = [spm_type(ctype) '_'];
    else
      prefix = job.prefix; 
    end
     
    if exist(V(1).fname,'file'), delete(V(1).fname); end
    if ndims(Y)==4
      %%
      N              = nifti;
      N.dat          = file_array(fullfile(pp,[prefix ff ee]),min([inf inf inf 3],size(Y)),[ctype spm_platform('bigend')],0,job.cvals,0);
      N.mat          = V(1).mat;
      N.mat0         = V(1).private.mat0;
      N.descrip      = V(1).descrip;
      create(N);    
      %Y(:,:,:,3) = Y(:,:,:,3) + Y(:,:,:,4);
      N.dat(:,:,:,:) = Y(:,:,:,:);
    else
      %%
      if isempty(job.prefix) && exist(fullfile(pp,[ff,ee]),'file'), delete(fullfile(pp,[ff,ee])); end
      V(1).fname    = fullfile(pp,[prefix ff ee]);
      V = rmfield(V,'private');
      spm_write_vol(V,Y);
    end
    spm_progress_bar('Set',si);
  end
  spm_progress_bar('Clear');
end