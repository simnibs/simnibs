function varargout = cat_vol_reduceRes(varargin)
% ______________________________________________________________________
% Reduction of the image resolution by merging of voxels with the name
% "lowresR#x#x#_*.nii with * for the original file name, and # for the 
% degree of reduction, with #=2 for the half resolutino.
%
% V = cat_vol_reduceRes(job)
% 
% job.data = cell of filename, use GUI if empty
% ______________________________________________________________________
% Robert Dahnke
% $Id: cat_vol_reduceRes.m 1036 2016-10-18 14:26:32Z dahnke $

  if nargin == 0 
      job.data = cellstr(spm_select([1 Inf],'image','select images to reduce'));
  else
      job = varargin{1};
  end

  if ~isfield(job,'data') || isempty(job.data)
     job.data = cellstr(spm_select([1 Inf],'image','select images to reduce'));
  else
     job.data = cellstr(job.data);
  end
  if isempty(job.data) || isempty(job.data{1})
    fprintf('No data!\n');  
    return;
  end
  
  def.verb   = 1; 
  def.lazy   = 1;
  def.prefix = 'lowres';
  job = cat_io_checkinopt(job,def);
  spm_clf('Interactive'); 
   
  V          = spm_vol(char(job.data));
  vx_vol     = sqrt(sum(V(1).mat(1:3,1:3).^2));
  if ~isfield(job,'res');
    job.res    = round(max(1,spm_input('Specify reduction factor(s)! ',...
                   '+1','i',round(min(vx_vol)*2.2 ./ vx_vol),[inf,3])));
  end
  if isempty(job.res)
    fprintf('No resolution!\n');  
    return;
  end
               
  Vr = repmat(V,1,size(job.res,1)); 
  Vr = rmfield(Vr,'private');
  
  spm_clf('Interactive'); 
  if job.verb, fprintf('cat_vol_reduceRes:\n'); end
  spm_progress_bar('Init',numel(job.data) .* size(job.res,1),'Resolution Reduction','Volumes Complete');
  for pi=1:numel(job.data)
    % load image
    Yi = single(spm_read_vols(V(pi))); 
    si = size(Yi); 
    
    %% processing
    for ri=1:size(job.res,1)
      % matrix size parameter
     
      ss = floor(si ./ job.res(ri,:));
      sx = ss .* job.res(ri,:) + (job.res(ri,:)-1); 
      if any(si<sx), Yi(sx(1),sx(2),sx(3))=0; end

      % filename
      [pp,ff,ee]       = spm_fileparts(job.data{pi});
      Vr(pi,ri).fname  = fullfile(pp,sprintf('%sR%dx%dx%d_%s%s',job.prefix,job.res(ri,:),ff,ee));
     
      if job.lazy && exist(Vr(pi,ri).fname,'file')
        Vr(pi,ri).dim   = ss;
        imat            = spm_imatrix(Vr(pi,ri).mat); 
        imat(1:3)       = imat(1:3) -  (job.res(ri,:)-1)/2;
        imat(7:9)       = imat(7:9) .* job.res(ri,:);
        Vr(pi,ri).mat   = spm_matrix(imat);
        if job.verb
          fprintf('Existing output %s\n',spm_file(Vr(pi,ri).fname,'link','spm_image(''%s'')'));
        end
        spm_progress_bar('Set',(pi-1) .* size(job.res,1) + ri); 
        continue 
      end
      
      %% image reduction 
      Yr = zeros(ss,'single'); 
      c  = 0;
      for x=1:job.res(ri,1)
        for y=1:job.res(ri,2)
          for z=1:job.res(ri,3)
            Yr = Yr + Yi( x:job.res(ri,1):job.res(ri,1)*ss(1) , ...
                          y:job.res(ri,2):job.res(ri,2)*ss(2) , ...
                          z:job.res(ri,3):job.res(ri,3)*ss(3) );
            c  = c + 1;
          end
        end
      end
      Yr = Yr / c; 

      %% write low resolution file
      Vr(pi,ri).dim   = size(Yr);
      imat            = spm_imatrix(Vr(pi,ri).mat); 
      imat(1:3)       = imat(1:3) -  (job.res(ri,:)-1)/2;
      imat(7:9)       = imat(7:9) .* job.res(ri,:);
      Vr(pi,ri).mat   = spm_matrix(imat);
      spm_write_vol(Vr(pi,ri),Yr);
      clear Y Yr imat;
    
      if job.verb
        fprintf('Output %s\n',spm_file(Vr(pi,ri).fname,'link','spm_image(''%s'')'));
      end
      spm_progress_bar('Set',(pi-1) .* size(job.res,1) + ri); 
    end
  end
  if job.verb, fprintf('cat_vol_reduceRes done.\n'); end
  spm_progress_bar('Clear');
  
  if nargout>0
    varargout{1} = Vr;
  end
end