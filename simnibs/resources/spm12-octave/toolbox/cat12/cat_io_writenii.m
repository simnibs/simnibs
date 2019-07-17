function varargout = cat_io_writenii(V,Y,folder,pre,desc,spmtype,range,writes,transform,YM,YMth)
% ______________________________________________________________________
% Write an image Y with the properties described by V with the datatype 
% spmtype for a specific range. Add the prefix pre and the description 
% desc to V. 
%
%   VO = cat_io_write_nii(Y,V[,pre,desc,spmtype,range,write,addpre,transform,YM,YMth])
%
%   Y       = input volume
%   V       = input volume structure
%   VO      = ouput volume structure
%   folder  = subfolder for writing data (default='')
%   pre     = prefix for filename (default='')
%   desc    = description that is added to the origin description (default='')
%   spmtype = spm image type (default given by the class of Y)
%   write   = [native warped modulated dartel]
%               native    0/1     (none/yes)
%               warped    0/1     (none/yes)
%               modulated 0/1/2/3 (none/affine+nonlinear/nonlinear only/both)
%               dartel    0/1/2/3 (none/rigid/affine/both)
%   transform = transformation data to write the image to warped, 
%               modulated, or dartel space (see cat_main)
%   YM      = mask for the final image (i.e. save thickness and ROIs)
%   YMth    = threshold for YM 
%
% Examples:
%   
%
% ______________________________________________________________________
% Robert Dahnke, Christian Gaser
% Structural Brain Mapping Group
% University Jena
%
% $Id: cat_io_writenii.m 1175 2017-09-06 14:45:25Z dahnke $
%
%#ok<*WNOFF,*WNON,*ASGLU>
  
  % file name
  if ~exist('pre','var'),  pre  = ''; end
  if ~exist('desc','var'), desc = ''; end

  % image type and convertations
  if ~exist('spmtype','var')
    switch class(Y)
      case 'logical',           spmtype = 'uint8';   
      case 'int8',              spmtype = 'int8';   
      case 'int16',             spmtype = 'int16'; 
      case 'int32',             spmtype = 'int32';  
      case {'uint8','char'},    spmtype = 'uint8'; 
      case 'uint16',            spmtype = 'uint16';
      case 'uint32',            spmtype = 'uint32'; 
      case {'single','double'}, spmtype = 'float32'; 
      otherwise
    end
  end  
  
  if ~exist('range','var'),  range  = [0 1]; end
  write = [1 0 0 0];
  if isstruct(writes)
    if isfield(writes,'native'),   write(1) = writes.native; end
    if isfield(writes,'warped'),   write(2) = writes.warped; end
    if isfield(writes,'mod'   ),   write(3) = writes.mod;    end
    if isfield(writes,'affine'),   write(4) = writes.affine; end
    if isfield(writes,'dartel'),   write(4) = writes.dartel; end
  elseif isnumeric(writes)
    if numel(writes)==3, write = [writes(1:2) 0 writes(3)]; else write = writes; end
  end
  if ~exist('YMth','var'), YMth = 0.5; end
  if exist('YM','var'),
    if all(size(Y)==size(YM))
      YM=single(YM);
    else
      error('MATLAB:cat_io_writenii:YM','Y and YM have different size');
    end
  end
  
  pp = fileparts(V.fname); 
  if ~exist('folder','var'), folder = ''; end
  if ~exist(fullfile(pp,folder),'dir'), mkdir(fullfile(pp,folder)); end
  
  
  % write native file
  % ____________________________________________________________________
  if write(1)==1
    fname = io_handle_pre(V.fname,pre,'',folder);
    if exist('transform','var') && isfield(transform,'native')
      if any(size(Y)~=transform.native.Vo.dim)
        nV = transform.native.Vi;
      else
        nV = transform.native.Vo;
      end
    else
      nV = V;
    end
    if exist(fname,'file'), delete(fname); end

    N         = nifti;
    N.dat     = file_array(fname,nV.dim(1:3),[spm_type(spmtype) ...
                  spm_platform('bigend')],range(1),range(2),0);
    N.mat     = nV.mat;

    % do not change mat0 - 20150612 
    % its a new image! 20150916
    warning off; 
    N.mat0    = nV.mat;
    warning on; 

    if isempty(V.descrip), N.descrip = desc; else  N.descrip = [desc ' < ' V.descrip]; end
    create(N);
    
    % final masking after transformation
    if exist('YM','var')
      switch spmtype
        case 'float32', N.dat(:,:,:) = double(Y) .* (smooth3(YM)>YMth); 
        otherwise,      N.dat(:,:,:) = max(range(1),min(range(1) + range(2)*double(intmax(spmtype)),double(Y)  .* (smooth3(YM)>YMth)));
      end
    else
      switch spmtype
        case 'float32', N.dat(:,:,:) = double(Y);
        otherwise,      N.dat(:,:,:) = max(range(1),min(range(1) + range(2)*double(intmax(spmtype)),double(Y)));
      end
    end

    Vn = spm_vol(fname); 
    % reduce to original native space if it was interpolated
    if exist('transform','var') && isfield(transform,'native') && any(size(Y)~=transform.native.Vo.dim)
      [pp,ff] = spm_fileparts(fname); 
      Vo = transform.native.Vo; 
      Vo.fname = fname; 
      Vo.dt    = Vn.dt; 
      Vo.pinfo = Vn.pinfo;
      if strcmp(ff(1:2),'p0')
        [Vn,Yn] = cat_vol_imcalc(Vn,Vo,'i1',struct('interp',5,'verb',0)); 
        % For Yp0 linear should be better, to avoid interpolation artifacts, but 
        % some single test showed that rounding produce better results.
        rf  = 100;
        Ynr = round(Yn*rf)/rf;
        YMR = false(size(Yn));
        for i=1:4, YMR = YMR | (Yn>(i-1/rf) & Yn<(i+1/rf)); end
        Yn(YMR)     = Ynr(YMR); clear YMR Ynr;
        delete(Vn.fname); % remove it, otherwise it will have the wrong filesize (correct readable, but still to big)
        Vn = spm_write_vol(Vn,double(Yn));
      else
        [Vn,Yn] = cat_vol_imcalc(Vn,Vo,'i1',struct('interp',6,'verb',0));
        delete(Vn.fname); % remove it, otherwise it will have the wrong filesize (correct readable, but still to big)
        Vn = spm_write_vol(Vn,double(Yn));
      end
    end
    
    if nargout>0, varargout{1}(1) = Vn; end
    if nargout>1, varargout{2}{1} = []; end
    clear N; 
  end

 
  % for masked images like thickness we need to fill undefined regions, 
  % to avoid the PVE of boundary voxel. 
  if any(write(2:end)) && exist('YM','var')
    [D,I] = cat_vbdist(single(Y)); Y(:)=Y(I(:)); clear D I; 
  end
  
  % deal with label maps 
  switch class(Y)
    case {'single','double'}
      labelmap = 0;
    case {'uint8','uint16'}
      if all(range == [0 1]); 
        labelmap = 1; 
        Y = single(Y); 
      else
        labelmap = 0;
      end
    otherwise
      labelmap = 0;
  end
  
  %% warped
  % ____________________________________________________________________
  % If we have a label map we have to correct the result, because spm_diffeo
  % and spm_field don't allow nearest neigbor deformation. Because the 
  % interpolated values of the boundaries cannot be simply rounded (it 
  % maybe generates another label), we need to replace this voxel by 
  % its nearest neighbor value.
  
  % interpolation to reduce artifacts if the resolution of the original 
  % image is similar or worse to that of the template resolution
  if write(2) || write(3)
    vx_vol   = sqrt(sum(transform.warped.M0(1:3,1:3).^2));  
    vx_volt  = sqrt(sum(transform.warped.M1(1:3,1:3).^2));  
    interpol = any(vx_vol>vx_volt) + any(vx_vol/2>vx_volt);
    if interpol
      if all(Y*1000 == round(Y*1000)) % label maps
        YI = interp3(Y,interpol,'nearest'); %yI = repmat(single(YI),1,1,1,3); 
      else
        YI = interp3(Y,interpol,'linear'); %yI = repmat(single(YI),1,1,1,3); 
      end
      %for i=1:3, yI(:,:,:,i) = interp3(transform.warped.yx(:,:,:,i),1,'linear'); end
    else
      YI = Y;
      %yI = transform.warped.y; 
    end
  end
  
  if write(2)
    pre2 = ['w'  pre]; desc2 = [desc '(warped)'];
    
    fname = io_handle_pre(V.fname,pre2,'',folder);
    if exist(fname,'file'), delete(fname); end
    if labelmap==0
      [wT,w]  = spm_diffeo('push',YI ,transform.warped.yx,transform.warped.odim(1:3));
      % divide by jacdet to get unmodulated data
      wT = wT./(w+0.001); 
    elseif labelmap==1
      % we can use modulated data throughout the following steps because the final maximum probability function
      % will be the same for modulated and unmodulated data
      wT = zeros([transform.warped.odim(1:3),max(YI(:))],'uint8'); 
      % interpolate each label separately
      for yi=1:max(YI(:)); 
        wTi  = spm_diffeo('push',single(YI==yi),transform.warped.yx,transform.warped.odim(1:3)); 
        wT(:,:,:,yi) = uint8(wTi*100); 
      end
      % use maximum probability function to get label again
      [wTmax,wT] = max(wT,[],4); 
    end  
    clear w;
    
    % final masking after transformation
    if exist('YM','var')
      [wTM,w] = spm_diffeo('push',YM,transform.warped.yx,transform.warped.odim(1:3));
      wTM = wTM./(w+0.001); 
      wTM = round(wTM*100)/100; 
      wT  = wT .* (smooth3(wTM)>YMth);
      clear w wTM;
    end
    
    N      = nifti;
    N.dat  = file_array(fname,transform.warped.odim, ...
              [spm_type(spmtype) spm_platform('bigend')], ...
              range(1),range(2),0);
    N.mat  = transform.warped.M1;
    N.mat0 = transform.warped.M1; % do not change mat0 - 20150612 - not changing, creating 20150916
    if isempty(V.descrip), N.descrip = desc; else  N.descrip = [desc2 ' < ' V.descrip]; end
    create(N);
    N.dat(:,:,:) = double(wT);
    clear N; 
    
    if nargout>0, varargout{1}(2) = spm_vol(fname); end
    if nargout>1, varargout{2}{2} = wT; end
  end

  
  
  %% modulated
  %  ___________________________________________________________________
  %{ 
  if write(3)
    %% 
    if write(3)==1 % SPM default with rigid transformation 
      pre3  = ['mw'   pre];
      desc3 = [desc '(Jac. sc. warped)'];
    elseif write(3)==2
      pre3 = ['m0w'  pre]; 
      desc3 = [desc '(Jac. sc. warped non-lin only)']; 
    end
    if (write(3)==1 && transform.warped.dartel==1) || (write(3)==2 && transform.warped.dartel==2)
    % for VBM Dartel, an affine registration was used and we now have to add the scaling factor

    elseif (write(3)==1 && transform.warped.dartel==2) || (write(3)==2 && transform.warped.dartel==1)
    % for VBM Shooting, a ridid registration was used and no further changes are required for standard SPM warped output

    end

    
    fname = io_handle_pre(V.fname,pre3,'',folder);
    if exist(fname,'file'), delete(fname); end
    
    [wT,wr] = spm_diffeo('push',YI,transform.warped.yx,transform.warped.odim(1:3)); 

    % final masking after transformation
    if exist('YM','var')
      wTM = spm_diffeo('push',YM,transform.warped.yx,transform.warped.odim(1:3)); 
      wT = wT .* (smooth3(wTM)>YMth);
    end
    
    % filtering of the jacobian determinant
    wrs = wr - 1; 
    spm_smooth(wrs,wrs,3/abs(transform.warped.M1(1))); wrs = wrs + 1;
    wT = spm_field(wr,wT ,[sqrt(sum(transform.warped.M1(1:3,1:3).^2)) 1e-6 1e-4 0  3 2]) .* wrs; 
    clear wrs;
  
    % create image
    N         = nifti;
    N.dat     = file_array(fname,transform.warped.odim,...
                  [spm_type(spmtype) spm_platform('bigend')], ...
                  range(1),range(2),0);
    N.mat     = transform.warped.M1;
    N.mat0    = transform.warped.M1; 
    if isempty(V.descrip), N.descrip = desc; else  N.descrip = [desc3 ' < ' V.descrip]; end
    create(N);       
    if write(3)==1
      N.dat(:,:,:) = double(wT)*abs(det(transform.warped.M0(1:3,1:3))/ ...
                     det(transform.warped.M1(1:3,1:3)));
    elseif write(3)==2
      N.dat(:,:,:) = double(wT)*abs(det(transform.warped.M2(1:3,1:3)));
    end
    clear N;
    
    if nargout>0, varargout{1}(3) = spm_vol(fname); end
    if nargout>1, varargout{2}{3} = wT*abs(det(transform.warped.M2(1:3,1:3))); end
  end
  %}
 
  %% old code ...
  if write(3)
    for wi=1:2  
      %% 
      if     (write(3)==1 || write(3)==3) && wi==1
        pre3 = ['mw'   pre]; desc3 = [desc '(Jac. sc. warped)'];
      elseif (write(3)==2 || write(3)==3) && wi==2
        pre3 = ['m0w'  pre]; desc3 = [desc '(Jac. sc. warped non-lin only)']; 
      else
        continue
      end

      fname = io_handle_pre(V.fname,pre3,'',folder);
      if exist(fname,'file'), delete(fname); end

      [wT,w]  = spm_diffeo('push',YI,transform.warped.yx,transform.warped.odim(1:3));

      % divide by jacdet to get unmodulated data
      wT = wT./(w+0.001); 

      % Modulation using spm_diffeo and push introduces aliasing artefacts,
      % thus we use the def2det function of the inverted deformations to obtain the old and 
      % in my view a more appropriate jacobian determinant 
      % The 2nd reason to use the old modulation is compatibility with cat_vol_defs.m
      Yy = spm_diffeo('invdef',transform.warped.yx,transform.warped.odim,eye(4),transform.warped.M0);
      w  = spm_diffeo('def2det',Yy)/det(transform.warped.M0(1:3,1:3));
      clear Yy

      % ensure that jacobian det is positive (no clue why some times the sign is switched)
      if mean(w(~isnan(w))) < 0, w = -w; end 
      w(:,:,[1 end]) = NaN; w(:,[1 end],:) = NaN; w([1 end],:,:) = NaN;
      wT = wT.*w;

      if interpol
        %[wTo,wo]  = spm_diffeo('push',Y,transform.warped.y,transform.warped.odim(1:3));
        %wT = wT * cat_stat_nansum(wTo(:))/cat_stat_nansum(wT(:)) * sum(Y(:))./sum(YI(:));
        wT = wT * cat_stat_nansum(Y(:))/cat_stat_nansum(wT(:)); 
      end
      clear w

      if exist('YM','var') % final masking after transformation
        if interpol
          YM = interp(YM,1,'linear'); 
        end
        wTM = spm_diffeo('push',YM,transform.warped.yx,transform.warped.odim(1:3)); 
        wT = wT .* (smooth3(wTM)>YMth);
        clear YM;
      end
     % clear yI

      % scale the jacobian determinant 
      if (write(3)==1 || write(3)==3) && wi==1
        wT = wT*abs(det(transform.warped.M0(1:3,1:3))/ ...
                    det(transform.warped.M1(1:3,1:3)));
      else
        wT = wT*abs(det(transform.warped.M2(1:3,1:3)));
      end

      N         = nifti;
      N.dat     = file_array(fname,transform.warped.odim, ...
                      [spm_type(spmtype) spm_platform('bigend')], ...
                      range(1),range(2),0);
      N.mat     = transform.warped.M1;
      N.mat0    = transform.warped.M1; % do not change mat0 - 20150612 - not changing, creating 20150916
      create(N);       
      if isempty(V.descrip), N.descrip = desc3; else  N.descrip = [desc3 ' < ' V.descrip]; end

      N.dat(:,:,:) = double(wT) ; % / 8^interpol; %
      clear N;

      if nargout>0, varargout{1}(3) = spm_vol(fname); end
      if nargout>1, varargout{2}{3} = wT; end
      clear wT
    end
  end
  
    
    
  %% write dartel files
  % ____________________________________________________________________
  if write(4)
    for wi=1:2
      if (write(4)==1 || write(4)==3) && wi==1 && isfield(transform,'rigid')
        transf=transform.rigid; 
        pre4=['r' pre]; post='_rigid'; desc4 = [desc ' (rigid)'];
      elseif (write(4)==2 || write(4)==3) && wi==2 && isfield(transform,'affine')
        transf=transform.affine;   
        pre4=['r' pre]; post='_affine'; desc4 = [desc ' (affine)'];
      else 
        continue
      end

     
      fname = io_handle_pre(V.fname,pre4,post,folder);
      if exist(fname,'file'), delete(fname); end
      VraT = struct('fname',fname,'dim',transf.odim,...
           'dt',   [spm_type(spmtype) spm_platform('bigend')],...
           'pinfo',[range(2) range(1)]','mat',transf.mat);%[1.0 0]'
      VraT = spm_create_vol(VraT);

      N  = nifti(VraT.fname);

      % do not change mat0 - 20150612
      % the mat0 contain the rigid transformation for the deformation tools!
      % get rid of the QFORM0 rounding warning
      warning off
      N.mat   = transf.mat; 
      N.mat0  = transf.mat0; 
      warning on

      %N.mat_intent  = 'Aligned';
      %N.mat0_intent = 'Aligned';
      if isempty(V.descrip), N.descrip = desc; else  N.descrip = [desc4 ' < ' V.descrip]; end
      create(N);

      for i=1:transf.odim(3),
        if labelmap
          tmp  = spm_slice_vol(double(Y) ,transf.M*spm_matrix([0 0 i]),transf.odim(1:2),0);
        else
          tmp  = spm_slice_vol(double(Y) ,transf.M*spm_matrix([0 0 i]),transf.odim(1:2),[1,NaN]);
        end
        if exist('YM','var')  % final masking after transformation
          if labelmap
            tmpM = spm_slice_vol(double(YM),transf.M*spm_matrix([0 0 i]),transf.odim(1:2),0);
          else
            tmpM = spm_slice_vol(double(YM),transf.M*spm_matrix([0 0 i]),transf.odim(1:2),[1,NaN]);
          end
          tmpM = smooth3(repmat(tmpM,1,1,3))>YMth;
          tmp  = tmp .* tmpM(:,:,2); clear tmpM; 
        end
        VraT = spm_write_plane(VraT,tmp,i);
      end

      if nargout>0, varargout{1}(4) = spm_vol(fname); end
      if nargout>1, varargout{2}{4} = []; end

      
    end
  end

  
end

function FO = io_handle_pre(F,pre,post,folder)
% Remove all known cat prefix types from a filename (and check if this file exist). 
  [pp,ff] = spm_fileparts(F); 

  % always use .nii as extension
  FO = fullfile(pp,folder,[pre ff post '.nii']);
end