function cat_io_3Dto4D(P,fname,avg)

%#ok<*ASGLU>

  if ~exist('P','var')
    P = spm_select(inf,'image','select image of the same space'); 
  end
  if isempty(P), return; end

  if ~exist('fname','var')
    fname = '4D';
  end
  
  if ~exist('avg','var')
    avg = 0;
  end
  
  V = spm_vol(P);
  P = cellstr(P);
  desc='';
  csv = {'ROIid','ROIname'};
  for i=1:numel(P)
    [pp,ff]=spm_fileparts(P{i}); 
    if any(V(1).dim(1:3) ~= V(i).dim(1:3)) %any(V(1).mat(:) ~= V(i).mat(:)) || 
      error('MATLAB:cat_io_3Dto4D:input_error','Bad resolution: %s %5d\n',i);
    end
    desc = sprintf('%s%s (%d),',ff,i);
    csv = [csv; {i,ff}]; %#ok<AGROW>
  end
  desc(end)='';
  cat_io_csv(fullfile(spm_fileparts(P{1}),[fname '.csv']),csv);
  
  if 1
    % real 4D-image
    N         = nifti;
    N.dat     = V(1).private.dat;
    N.dat.fname  = fullfile(spm_fileparts(P{1}),[fname '.nii']);
    N.dat.dim(4) = numel(P);
    N.mat     = V(1).mat;
    N.mat0    = V(1).private.mat0;
    N.descrip = desc;
    create(N);       

    for i=1:numel(P)
      Y = spm_read_vols(V(i));
      N.dat(:,:,:,i) = Y;
    end
  end
  
  if avg
    N         = nifti;
    N.dat     = V(1).private.dat;
    N.dat.fname = fullfile(spm_fileparts(P{1}),'A4D.nii');
    N.mat     = V(1).mat;
    N.mat0    = V(1).private.mat0;
    N.descrip = desc;
    create(N);       
    Y = spm_read_vols(spm_vol(fullfile(spm_fileparts(P{1}),[fname '.nii']))); 
    [maxx,N.dat(:,:,:)] = nanmax(Y,[],4);
  end
end