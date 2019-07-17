function cat_run_job(job,tpm,subj)
% run CAT 
% ______________________________________________________________________
%
% Initialization function of the CAT preprocessing. 
%  * creation of the subfolder structure (if active)
%  * check of image resolution (avoid scans with very low resolution)
%  * noise correction (ISARNLM)
%  * interpolation 
%  * affine preprocessing (APP)
%    >> cat_run_job_APP_SPM
%    >> cat_run_job_APP_init
%    >> cat_run_job_APP_final
%  * affine registration
%  * initial SPM preprocessing
%
%   cat_run_job(job,tpm,subj)
% 
%   job  .. SPM job structure with main parameter
%   tpm  .. tissue probability map (hdr structure)
%   subj .. file name
% ______________________________________________________________________
% Christian Gaser
% $Id: cat_run_job.m 1245 2017-12-12 10:35:00Z dahnke $

%#ok<*WNOFF,*WNON>

    % if there is a breakpoint in this file set debug=1 and do not clear temporary variables 
    dbs   = dbstatus; debug = 0; for dbsi=1:numel(dbs), if strcmp(dbs(dbsi).name,mfilename); debug = 1; break; end; end

    global cat_err_res; % for CAT error report

    stime = clock;

    
    % create subfolders if not exist
    pth = spm_fileparts(job.channel(1).vols{subj}); 
    if job.extopts.subfolders
      folders = {'mri','report'};
      for i=1:numel(folders)
        if ~exist(fullfile(pth,folders{i}),'dir')
          mkdir(fullfile(pth,folders{i}));
        end
      end
      if ~exist(fullfile(pth,'surf'),'dir') && job.output.surface
        mkdir(fullfile(pth,'surf'));
      end
      if ~exist(fullfile(pth,'label'),'dir') && job.output.ROI
        mkdir(fullfile(pth,'label'));
      end
      
      mrifolder    = 'mri';
      reportfolder = 'report';
    else
      mrifolder    = '';
      reportfolder = '';
    end
  
    
    % create subject-wise diagy file with the command-line output
    [pp,ff,ee,ex] = spm_fileparts(job.data{subj}); 
    catlog = fullfile(pth,reportfolder,['catlog_' ff '.txt']);
    if exist(catlog,'file'), delete(catlog); end % write every time a new file, turn this of to have an additional log file
    diary(catlog); 
    
    
    % print current CAT release number and subject file
    [n,r] = cat_version;
    str  = sprintf('CAT12 r%s: %d/%d',r,subj,numel(job.channel(1).vols));
    str2 = spm_str_manip(job.channel(1).vols{subj}(1:end-2),['a' num2str(70 - length(str))]);
    cat_io_cprintf([0.2 0.2 0.8],'\n%s\n%s: %s%s\n%s\n',...
          repmat('-',1,72),str,...
          repmat(' ',1,70 - length(str) - length(str2)),str2,...
          repmat('-',1,72));
    clear r str str2


    
    %  -----------------------------------------------------------------
    %  separation of full CAT preprocessing and SPM segmentation
    %  preprocessing (running DARTEL and PBT with SPM segmentation)
    %  -----------------------------------------------------------------
    [pp,ff,ee,ex] = spm_fileparts(job.data{subj}); 
    if exist(fullfile(pp,['c1' ff(3:end) ee]),'file') && ...
       exist(fullfile(pp,['c2' ff(3:end) ee]),'file') && ...
       exist(fullfile(pp,['c3' ff(3:end) ee]),'file') && ...
       exist(fullfile(pp,[ff(3:end) '_seg8.mat']),'file');
       
        job.data{subj}          = fullfile(pp,[ff ee]); 
        job.channel.vols{subj}  = fullfile(pp,[ff ee]); 

        % prepare SPM preprocessing structure 
        images = job.channel(1).vols{subj};
        for n=2:numel(job.channel)
          images = char(images,job.channel(n).vols{subj});
        end

        obj.image    = spm_vol(images);
        spm_check_orientations(obj.image);

        obj.fwhm     = job.opts.fwhm;
        obj.biasreg  = cat(1,job.opts.biasreg);
        obj.biasfwhm = cat(1,job.opts.biasfwhm);
        obj.tpm      = tpm;
        obj.lkp      = [];
        if all(isfinite(cat(1,job.tissue.ngaus))),
            for k=1:numel(job.tissue),
                obj.lkp = [obj.lkp ones(1,job.tissue(k).ngaus)*k];
            end;
        end
        obj.reg      = job.opts.warpreg;
        obj.samp     = job.opts.samp;              

        cfname  = fullfile(pp,[ff ee]);
        ofname  = fullfile(pp,[ff(3:end) ee]); 
        nfname  = fullfile(pp,mrifolder,['n' ff '.nii']); 
        copyfile(ofname,nfname); 

        res = load(fullfile(pp,[ff(3:end) '_seg8.mat']));
        job.channel(1).vols{subj}  = [nfname ex];
        job.channel(1).vols0{subj} = [ofname ex];
        res.image  = spm_vol([nfname ex]);
        res.image0 = spm_vol([ofname ex]);
        res.imagec = spm_vol([cfname ex]);
        res.spmpp  = 1; 
    else

        %  -----------------------------------------------------------------
        %  check resolution properties
        %  -----------------------------------------------------------------
        %  There were some images that should not be processed. So we have  
        %  to check for high slice thickness and low resolution.
        %  -----------------------------------------------------------------
        for n=1:numel(job.channel) 
          V = spm_vol(job.channel(n).vols{subj});
          vx_vol = sqrt(sum(V.mat(1:3,1:3).^2));

          if any(vx_vol>5)  % too thin slices
            error('CAT:cat_main:TooLowResolution', sprintf(...
                 ['Voxel resolution has to be better than 5 mm in any dimension \n' ...
                  'for reliable CAT preprocessing! \n' ...
                  'This image has a resolution %0.2fx%0.2fx%0.2f mm%s. '], ... 
                    vx_vol,char(179))); %#ok<SPERR>
          end
          if prod(vx_vol)>27  % too small voxel volume (smaller than 3x3x3 mm3)
            error('CAT:cat_main:TooHighVoxelVolume', ...
                 ['Voxel volume has to be smaller than 10 mm%s (around 3x3x3 mm%s) to \n' ...
                  'allow a reliable CAT preprocessing! \n' ...
                  'This image has a voxel volume of %0.2f mm%s. '], ...
                  char(179),char(179),prod(vx_vol),char(179));
          end
          if max(vx_vol)/min(vx_vol)>8 % isotropy 
            error('CAT:cat_main:TooStrongIsotropy', sprintf(...
                 ['Voxel isotropy (max(vx_size)/min(vx_size)) has to be smaller than 8 to \n' ...
                  'allow a reliable CAT preprocessing! \n' ...
                  'This image has a resolution %0.2fx%0.2fx%0.2f mm%s and a isotropy of %0.2f. '], ...
                  vx_vol,char(179),max(vx_vol)/min(vx_vol))); %#ok<SPERR>
          end
        end


        % save original file name 
        for n=1:numel(job.channel) 
          job.channel(n).vols0{subj} = job.channel(n).vols{subj};
        end

        % allways create the n*.nii image because of the real masking of the
        % T1 data for spm_preproc8 that include rewriting the image!
        for n=1:numel(job.channel) 
            [pp,ff,ee] = spm_fileparts(job.channel(n).vols{subj}); 
            ofname  = fullfile(pp,[ff ee]); 
            nfname  = fullfile(pp,mrifolder,['n' ff '.nii']); 
            if strcmp(ee,'.nii')
              copyfile(ofname,nfname); 
            elseif strcmp(ee,'.img')
              V = spm_vol(job.channel(n).vols{subj});
              Y = spm_read_vols(V);
              V.fname = nfname; 
              spm_write_vol(V,Y);
              clear Y; 
            end
            job.channel(n).vols{subj} = nfname;

            
            
            %% skull-stripping detection
            %  ------------------------------------------------------------
            %  Detect skull-stripping or defaceing because it strongly 
            %  affects SPM segmenation that expect gaussian distribution! 
            %  If a brain mask was used than we expect 
            %   - many zeros (50% for small background - 80-90% for large backgrounds)
            %   - a brain like volume (below 2500 cm3)
            %   - only on object (the masked regions)
            %   - only on background (not in very case?)
            %   - less variance of thissue intensity (only 3 brain classes)
            %  ------------------------------------------------------------
            VF    = spm_vol(nfname); 
            YF    = spm_read_vols(VF); 
            Oth   = cat_stat_nanmean(YF(YF(:)~=0 & YF(:)>cat_stat_nanmean(YF(:)))); 
            F0vol = cat_stat_nansum(YF(:)~=0) * prod(vx_vol) / 1000; 
            F0std = cat_stat_nanstd(YF(YF(:)>0.5*Oth & YF(:)>0)/Oth); 
            YFC = YF~=0; 
            if sum(YFC(:)>0)<numel(YFC)*0.9 && sum(YFC(:)>0)>numel(YFC)*0.1  % if there is a meanful background
              YFC = ~cat_vol_morph(YF~=0,'lc',1);                            % close noisy background
            end
            [YL,numo] = spm_bwlabel(double(YF~=0),26);  clear YL;            % number of objects
            [YL,numi] = spm_bwlabel(double(YFC==0),26); clear YL;            % number of background regions 
            ppe.affreg.skullstrippedpara = [sum(YF(:)==0)/numel(YF) numo numi F0vol F0std]; 
            ppe.affreg.skullstripped = ...
              ppe.affreg.skullstrippedpara(1)>0.5 && ...                     % many zeros
              ppe.affreg.skullstrippedpara(2)<5  && ...                      % only few object
              ppe.affreg.skullstrippedpara(3)<10 && ...                      % only few background regions 
              F0vol<2500 && F0std<0.5;                                       % many zeros  and  not to big
            ppe.affreg.skullstripped = ppe.affreg.skullstripped || ...
              sum([ppe.affreg.skullstrippedpara(1)>0.8 F0vol<1500 F0std<0.4])>1; % or 2 extrem values
            if ~debug, clear YFC F0vol F0std numo numi; end 
    
            
            
            %% APP bias correction (APP1 and APP2)
            %  ------------------------------------------------------------
            %  Bias correction is essential for stable affine registration. 
            %  SPM further required Gaussian distributed data that is 
            %  achieved by Smoothing in high resolution data and by 
            %  additional noise in regions with many zeros typical in 
            %  skull-stripped or defaced data. 
            %  ------------------------------------------------------------
            if (job.extopts.APP==1 || job.extopts.APP==2) && ~ppe.affreg.skullstripped
              stime = cat_io_cmd('APPs bias correction');
              if job.extopts.verb>1, fprintf('\n'); end
              stime2 = cat_io_cmd('  Preparation','g5','',job.extopts.verb-1);
              if debug, copyfile(ofname,nfname); end
              
              % SPM segmentation parameter
              preproc.channel.vols{1,1}   = nfname; 
              preproc.channel.biasreg     = min(0.01,max(0.0001,job.opts.biasreg));
              preproc.channel.biasfwhm    = min(90,max(30,job.opts.biasfwhm/2));
              preproc.channel.write       = [0 1];
              for ti=1:6
                preproc.tissue(ti).tpm    = {[tpm.V(ti).fname ',' num2str(ti)]};
                preproc.tissue(ti).ngaus  = job.opts.ngaus(ti);
                preproc.tissue(ti).native = [0 0];                           % native, dartel
                preproc.tissue(ti).warped = [0 0];                           % unmod, mod
              end
              preproc.warp.mrf            = 1;                               % 
              preproc.warp.cleanup        = 1;                               % this is faster and give better results!
              preproc.warp.reg            = [0 0.001 0.5 0.05 0.2];
              preproc.warp.affreg         = job.opts.affreg; 
              preproc.warp.write          = [0 0];
              %preproc.warp.fwhm           = 0;
              %preproc.warp.samp           = min(9,max(2,job.opts.samp*2));
              
              % add noise in zerosed regions (skull-stripping / defacing)
              VF = spm_vol(nfname);
              YF = single(spm_read_vols(VF));
              % some average object intensity 
              Tthn = cat_stat_nanmean(YF(YF(:)>cat_stat_nanmean(YF(YF(:)~=0)) & YF(:)~=0));
              % limitation is required for division data (yv98_05mm_corrected_
              YF = min( Tthn*10, YF); 
              % smoothing for Gaussian distribution
              YF(isnan(YF(:)))=0;
              YFs = YF+0; spm_smooth(YFs,YFs,0.5./vx_vol);
              YM = abs( (YF-YFs)./max(eps,YFs));
              YM = min(1,smooth3(YM ./ min(0.2,max(YM(:)))));
              YF = YFs.*YM + YF.*(1-YM);
              Y0 = YF==0; 
              % add noise for Gaussian distribution in the background
              if ppe.affreg.skullstripped
                YF = YF + (Y0) .* (0.05*Tthn).*rand(size(YF));
              else
                YF(cat_vol_morph(Y0,'o'))=nan; 
              end
              % force floating point
              VF.dt(1) = 16; 
              spm_write_vol(VF,YF);  
              clear VF YF Tthn; 

            
              %% try SPM preprocessing
              %  ----------------------------------------------------------
              %  * if SPM failed we go on without bias correction 
              %  * further interations on different resolution are usefull
              %    (the first run correct the worst problems and should 
              %     a correct registration of the next runs)
              %  * best results for strong corrections (fwhm 30 to 45 mm)
              %  * SPM preprocessing is very fast (and gives you results in 
              %    the template resolution?) but writing results in the 
              %    original resolution is very slow, so using 2 interations
              %    is maybe optimal
              %  * much slower for biasfwhm<35 mm and samp<4.5 mm
              %  * the final operation is much slower (3 times) because it
              %    required the estimation of the deformation field to
              %    write data in the orignal space
              %  * writing further outputs on disk will also cost more time
              %    >> further optimization is possible by avoiding
              %       temporary disk output 
              %  ----------------------------------------------------------
              if job.opts.biasfwhm<45                               % strong  (3 iterations)
                fwhmx  = 2:-1/4:3/4; 
                sampx  = 2:-0.5/(numel(fwhmx)-1):1.5; %repmat(2,numel(fwhmx),1); 
              elseif job.opts.biasfwhm>=45 && job.opts.biasfwhm<75      % medium (2 iterations) 
                fwhmx  = 2:-1/4:3/4; 
                sampx  = 2:-0.5/(numel(fwhmx)-1):1.5; %repmat(2,numel(fwhmx),1); 
              elseif job.opts.biasfwhm>=75                              % light (1 iteration)
                factor = 1/2;
              end
              %
              %fwhmx  = 2:-1/4:3/4; 
              %sampx  = 2:-0.5/(numel(fwhmx)-1):1.5; %repmat(2,numel(fwhmx),1); 

                
                
              spmp0    = debug;  % job.extopts.verb>1; % for debugging: 0 - remove all APP data, 1 - save Ym, 2 - save Ym and Yp0
              optimize = 0;      % use low resolution (in mm) input image to increase speed >> limited speed advante :/
                                 % deformation requires a lot of time ... maybe its enough to use it only in case of ultra highres data
              if any(vx_vol<1.5), optimize = 1.5; end            
              
              preproc.Yclsout = false(1,6); 
              if debug, copyfile(ofname,nfname); end
              if optimize>0 
                %% lower resolution to improve SPM processing time and use smoothing for denoising
                Vi        = spm_vol(nfname); 
                Vi        = rmfield(Vi,'private'); Vn = Vi; 
                imat      = spm_imatrix(Vi.mat); 
                Vi.dim    = round(Vi.dim .* vx_vol./repmat(optimize,1,3));
                imat(7:9) = repmat(optimize,1,3) .* sign(imat(7:9));
                Vi.mat    = spm_matrix(imat);

                spm_smooth(nfname,nfname,repmat(0.75,1,3)); % denoising
                [Vi,Yi] = cat_vol_imcalc(Vn,Vi,'i1',struct('interp',1,'verb',0));
                delete(nfname); spm_write_vol(Vi,Yi); 
                %if ~isinf(job.opts.biasstr), clear Yi; end
              end
              
              %%
              bias=0; 
              for ix=1:numel(sampx)
                %% parameter update
                preproc.warp.samp         = min(9 ,max(1 ,job.opts.samp     * sampx(ix))); 
                preproc.channel.biasfwhm  = min(120,max(30,job.opts.biasfwhm * fwhmx(ix)));
                % preproc.warp.reg          = [0 0 0 0 0];
      
                stime2 = cat_io_cmd(sprintf('  SPM bias correction (samp: %0.2f mm, fwhm: %3.0f mm)',preproc.warp.samp,...
                  preproc.channel.biasfwhm),'g5','',job.extopts.verb-1,stime2);
                
                %try
                  % SPM bias correction
                  warning off; 
                  if ix==numel(sampx) || bias(1) >= preproc.channel.biasfwhm 
                    if job.extopts.APP==2 || spmp0>1, preproc.Yclsout = true(1,6); end % 
                    vout   = cat_spm_preproc_run(preproc,'run');
                    spmmat = open(fullfile(pp,mrifolder,['n' ff '_seg8.mat'])); 
                    AffineSPM = spmmat.Affine; 
                  else
                    cat_spm_preproc_run(preproc,'run');
                  end
                  warning on; 

                  Pmn    = fullfile(pp,mrifolder,['mn' ff ee]); 
                  Pmn_ri = fullfile(pp,mrifolder,['mn' ff '_r' num2str(ix) ee]);
                  Pmn_r0 = fullfile(pp,mrifolder,['mn' ff '_r0' ee]);
                  
                  % estimate bias strength based on the applied corrections
                  % of the initial correction .. 
                  % in case of updates, local biasfield strenght is maybe
                  % better (only useful if strong changes are allowed)
                  if ix==1 
                    Vn   = spm_vol(Pmn); Vn = rmfield(Vn,'private'); Yn = spm_read_vols(Vn);
                    bias(ix) = (1/cat_stat_nanstd(Yn(:)./Yi(:))) * 4; 
                    fprintf('bias=%5.0f mm ',bias(ix)); 
                    bias(ix) = max(30,min(120,round(bias(ix) / 15) * 15));
                    if ~debug, clear Yn; end
                  %elseif ix>1 && isinf(job.opts.biasstr)
                    % nothing to do
                  else 
                    bias(ix) = 0; 
                  end
                  
                  % backup the bias corrected image
                  if spmp0>0, copyfile(Pmn,Pmn_ri); end
                  
                  % backup for mixing
                  if (ix==2 && numel(sampx)>2) || (ix==1 && numel(sampx)<=2) || bias(1) > preproc.channel.biasfwhm 
                    copyfile(fullfile(pp,mrifolder,['mn' ff ee]),Pmn_r0); 
                  end
                  copyfile(fullfile(pp,mrifolder,['mn' ff ee]),nfname);
                  
                  % write segmentation 
                  if spmp0>1 && exist('vout','var') && isfield(vout,'Ycls') 
                    VF  = spm_vol(nfname); VF.fname = fullfile(pp,mrifolder,['p0n' ff '.nii']); VF.dt(1) = 2; 
                    Yp0 = single(vout.Ycls{1})/255*2 + single(vout.Ycls{2})/255*3 + single(vout.Ycls{3})/255;
                    spm_write_vol(VF,Yp0);  
                  elseif spmp0==0 && ix==numel(sampx) % remove spm mat file
                    delete(fullfile(pp,mrifolder,['n' ff '_seg8.mat']));
                  end
                  
                  
                  
                  %% combine low and high frequency filted images 
                  if ix==numel(sampx) || (bias(1) > preproc.channel.biasfwhm && ix>1)
                    %%
                    stime2 = cat_io_cmd('  Postprocessing','g5','',job.extopts.verb-1,stime2);
                    if optimize
                      % update Ym
                      Vn  = spm_vol(nfname); Vn = rmfield(Vn,'private');
                      Vo  = spm_vol(ofname); Vo = rmfield(Vo,'private'); 
                      [Vmx,vout.Ym] = cat_vol_imcalc([Vo,Vn],Vo,'i2',struct('interp',1,'verb',0));
                      vout.Ym = single(vout.Ym); 

                      %% remap Ym0
                      Vn = spm_vol(Pmn_r0); Vn = rmfield(Vn,'private');
                      Vo = spm_vol(ofname); Vo = rmfield(Vo,'private'); 
                      [Vmx,Ym0] = cat_vol_imcalc([Vo,Vn],Vo,'i2',struct('interp',1,'verb',0));
                      Ym0 = single(Ym0); 

                      %% update Ycls
                      if isfield(vout,'Ycls')
                        for i=1:6
                          Vn  = spm_vol(nfname); Vn = rmfield(Vn,'private');
                          Vo  = spm_vol(ofname); Vo = rmfield(Vo,'private'); 
                          Vn.pinfo = repmat([1;0],1,size(vout.Ycls{i},3));
                          Vo.pinfo = repmat([1;0],1,Vo.dim(3));
                          Vn.dt    = [2 0];
                          Vo.dt    = [2 0]; 
                          Vo.dat   = zeros(Vo.dim(1:3),'uint8');
                          if ~isempty(vout.Ycls{i})
                            Vn.dat   = vout.Ycls{i}; 
                            [Vmx,vout.Ycls{i}]  = cat_vol_imcalc([Vo,Vn],Vo,'i2',struct('interp',1,'verb',0));
                            vout.Ycls{i} = cat_vol_ctype(vout.Ycls{i}); 
                          end
                        end
                      end
                    else
                      Ym0 = spm_read_vols(spm_vol(Pmn_r0)); 
                    end
                    if debug, stime2 = cat_io_cmd(' ','g5','',job.extopts.verb-1,stime2); end
                    
                    if ~debug && exist(Pmn_r0,'file'), delete(Pmn_r0); end
                    if ~debug && exist(Pmn_r0,'file'), delete(Pmn_r0); end
                    
                    %% mixing
                    %  creation of segmenation takes a lot of time because
                    %  of the deformations. So it is much faster to load a
                    %  rought brain mask. 
                    if isfield(vout,'Ycls')
                      Yp0 = single(vout.Ycls{1})/255*2 + single(vout.Ycls{2})/255*3 + single(vout.Ycls{3})/255;
                      YM2 = cat_vol_smooth3X(cat_vol_smooth3X(Yp0>0,16/mean(vx_vol))>0.95,10/mean(vx_vol)); 
                      if ~debug, clear Yp0; end
                    else
                      Pb  = char(job.extopts.brainmask);
                      Pbt = fullfile(pp,mrifolder,['brainmask_' ff '.nii']);
                      VF  = spm_vol(ofname); 
                      VFa = VF; %if job.extopts.APP~=5, VFa.mat = Affine * VF.mat; end
                      [Vmsk,Yb] = cat_vol_imcalc([VFa,spm_vol(Pb)],Pbt,'i2',struct('interp',3,'verb',0)); 
                      Ybb = cat_vol_smooth3X(Yb>0.5,8/mean(vx_vol)); Ybb = Ybb./max(Ybb(:)); 
                      YM2 = cat_vol_smooth3X(Ybb>0.95,8/mean(vx_vol)); YM2 = YM2./max(Ybb(:)); 
                      if ~debug, clear Pb Pbt VFa Vmsk Yb Ybb; end
                    end
                    %% combine the low (~60 mm, Ym0) and high frequency correction (~30 mm, vout.Ym) 
                    Yo  = spm_read_vols(spm_vol(ofname)); 
                    Yw  = vout.Ym.*(1-YM2) + (YM2).*Ym0;
                    Yw  = Yo./Yw .* (Yo~=0 & Yw~=0); 
                    % correct undefined voxel and assure smoothness of the bias field 
                    Yw  = cat_vol_approx(max(0,min(2,Yw)),'',vx_vol,2,struct('lfO',2));
                    vout.Ym = Yo ./ Yw; 
                    if ~debug, clear Yw Yo; end
                    Vm = spm_vol(ofname); Vm.fname = nfname; Vm.dt(1) = 16;  
                    spm_write_vol(Vm,vout.Ym); 
                    
                    % backup the bias corrected image
                    if spmp0>0 
                      copyfile(nfname,fullfile(pp,mrifolder,['mn' ff '_r' num2str(ix)+1 ee])); 
                    end
                    if exist(Pmn_r0,'file'), delete(Pmn_r0); end
                    break
                  else
                    movefile(fullfile(pp,mrifolder,['mn' ff ee]),nfname); 
                  end
                try
                catch
                  fprintf('\b\b\b\b\b\b\b\b\b(failed) ');   
                  if exist(fullfile(pp,mrifolder,['mn' ff ee]),'file')
                    delete(fullfile(pp,mrifolder,['mn' ff ee]));
                  end
                end
                if 0
                  %% just debugging
                  Yp0  = single(vout.Ycls{1})/255*2 + single(vout.Ycls{2})/255*3 + single(vout.Ycls{3})/255;
                  Tthx = max( [ median( vout.Ym(vout.Ycls{1}(:)>128)) , median( vout.Ym(vout.Ycls{2}(:)>128)) ]);
                  ds('d2','a',vx_vol,Ym0/Tthx,Yp0/3,vout.Ym/Tthx,vout.Ym/Tthx,120)
                  %%
                  Tthx = mean( vout.Ym( YM2(:)>0.2 & vout.Ym(:)>mean(vout.Ym(YM2(:)>0.2))) );
                  ds('d2','',vx_vol,Ym0/Tthx,YM2,vout.Ym/Tthx,Ym/Tthx,80)
                end
              end
              if debug, cat_io_cmd('','g5','',job.extopts.verb-1,stime); end
              Pmn = fullfile(pp,mrifolder,['mn' ff ee]); 
              if exist(Pmn,'file'), delete(Pmn); end
              
              %% try APPs preprocessing
              %  ----------------------------------------------------------
              %  SPM may failed in correction of high intensive gyri that 
              %  should be fixed by the APPs correction that is similar to
              %  the APPinit and APPfinal routines but used SPM as input 
              %  segmentation.
              %  However this can work very well in general, problems can 
              %  occure in subcortical structures and the cerebellum, 
              %  especially in special protocols. 
              %  It is also possible to skip this function if the estimated
              %  bias field is very small, but I am not sure if the is useful.
              %  ----------------------------------------------------------
              if exist('vout','var') && job.extopts.APP==2 && isfield(vout,'Ycls')
                stime2 = cat_io_cmd('  APP bias correction','g5','',job.extopts.verb-1,stime2); 
                
                try
                  if isinf(job.opts.biasstr), catbias = 1 - (bias(1)-30)/60; else catbias = job.opts.biasstr; end
                  Ym = cat_run_job_APP_SPM(ofname,vout,vx_vol * (2 - strcmp(job.extopts.species,'human')),...
                    job.extopts.verb-1,catbias); 
                  spm_write_vol(spm_vol(nfname),Ym);  
                  if ~debug, clear vout Ym0 YM2 Ym; end  
                  if spmp0
                    copyfile(nfname,fullfile(pp,mrifolder,['mn' ff '_r' num2str(ix)+2 ee])); 
                  end
                catch
                  fprintf('failed\n');
                end
              end
              fprintf('%5.0fs\n',etime(clock,stime));  
              
              %if spmp0, return; else if ~debug, clear spmp0; end; end
            end
            
            
            
            %% noise correction
            %  ------------------------------------------------------------
            if job.extopts.NCstr~=0
              if job.extopts.NCstr==2 || job.extopts.NCstr==3
                if job.extopts.NCstr==2, NCstr=-inf; else NCstr=1; end 
                stime = cat_io_cmd(sprintf('ISARNLM denoising (NCstr=%0.2f)',NCstr));
                if job.extopts.verb>1, fprintf('\n'); end
                cat_vol_isarnlm(struct('data',nfname,'verb',(job.extopts.verb>1)*2,'prefix','','NCstr',NCstr)); 
                if job.extopts.verb>1, cat_io_cmd(' ','',''); end
              else
                stime = cat_io_cmd(sprintf('SANLM denoising (NCstr=%0.2f)',job.extopts.NCstr));
                cat_vol_sanlm(struct('data',nfname,'verb',0,'prefix','','NCstr',job.extopts.NCstr)); 
              end
              fprintf('%5.0fs\n',etime(clock,stime));   
            end
        end



        %% Interpolation
        %  -----------------------------------------------------------------
        %  The interpolation can help to reduce problems for morphological
        %  operations for low resolutions and strong isotropic images. 
        %  For interpolation from 2 to 1 mm the segmentation shows more 
        %  anatomical details.
        for n=1:numel(job.channel) 

          % prepare header of resampled volume
          Vi        = spm_vol(job.channel(n).vols{subj}); 
          vx_vol    = sqrt(sum(Vi.mat(1:3,1:3).^2));
          vx_vol    = round(vx_vol*10^2)/10^2; % avoid small differences 

          % we have to look for the name of the field due to the GUI job struct generation! 
          restype   = char(fieldnames(job.extopts.restypes));
          switch restype
            case 'native'
              vx_voli  = vx_vol;
            case 'fixed', 
              vx_voli  = min(vx_vol ,job.extopts.restypes.(restype)(1) ./ ...
                         ((vx_vol > (job.extopts.restypes.(restype)(1)+job.extopts.restypes.(restype)(2)))+eps));
              vx_voli  = max(vx_voli,job.extopts.restypes.(restype)(1) .* ...
                         ( vx_vol < (job.extopts.restypes.(restype)(1)-job.extopts.restypes.(restype)(2))));
            case 'best'
              best_vx  = max( min(vx_vol) ,job.extopts.restypes.(restype)(1)); 
              vx_voli  = min(vx_vol ,best_vx ./ ((vx_vol > (best_vx + job.extopts.restypes.(restype)(2)))+eps));
            otherwise 
              error('cat_run_job:restype','Unknown resolution type ''%s''. Choose between ''fixed'',''native'', and ''best''.',restype)
          end


          % interpolation 
          if any( (vx_vol ~= vx_voli) )  
            stime = cat_io_cmd(sprintf('Internal resampling (%4.2fx%4.2fx%4.2fmm > %4.2fx%4.2fx%4.2fmm)',vx_vol,vx_voli));
           
            Vi        = rmfield(Vi,'private');
            imat      = spm_imatrix(Vi.mat); 
            Vi.dim    = round(Vi.dim .* vx_vol./vx_voli);
            imat(7:9) = vx_voli .* sign(imat(7:9));
            Vi.mat    = spm_matrix(imat);

            Vn = spm_vol(job.channel(n).vols{subj}); 
            Vn = rmfield(Vn,'private'); 
            cat_vol_imcalc(Vn,Vi,'i1',struct('interp',2,'verb',0));
            vx_vol = vx_voli;
          
            fprintf('%5.0fs\n',etime(clock,stime));     
          else
            vx_vol = sqrt(sum(Vi.mat(1:3,1:3).^2));
          end
          clear Vi Vn;
        end


        %  prepare SPM preprocessing structure 
        images = job.channel(1).vols{subj};
        for n=2:numel(job.channel)
            images = char(images,job.channel(n).vols{subj});
        end
        obj.image    = spm_vol(images);
        obj.fwhm     = job.opts.fwhm;
        obj.biasreg  = job.opts.biasreg;
        obj.biasfwhm = job.opts.biasfwhm; 
        obj.tpm      = tpm;
        obj.lkp      = [];
        if all(isfinite(cat(1,job.tissue.ngaus))),
            for k=1:numel(job.tissue),
                obj.lkp = [obj.lkp ones(1,job.tissue(k).ngaus)*k];
            end;
        end
        obj.reg      = job.opts.warpreg;
        obj.samp     = job.opts.samp;              
        spm_check_orientations(obj.image);



      
        if 1 || ~exist('Affine','var')
          %% Initial affine registration.
          %  -----------------------------------------------------------------
          %  APP option with subparameter
          %  Skull-stripping is helpful for correcting affine registration of neonates and other species. 
          %  Bias correction is important for the affine registration.
          %  However, the first registation can fail and further control is required   
          if exist('AffineSPM','var'), Affine = AffineSPM; else Affine  = eye(4); end
          [pp,ff] = spm_fileparts(job.channel(1).vols{subj});
          Pbt     = fullfile(pp,mrifolder,['brainmask_' ff '.nii']);
          Pb      = char(job.extopts.brainmask);
          Pt1     = char(job.extopts.T1);


          % Estimate the COM and initial affine matrix that is used in case 
          % of failed initial affine registration.
          VF    = spm_vol(nfname); 
          YF    = spm_read_vols(VF); 
          xsum  = shiftdim(sum(sum(YF/Oth>0.5,2),3),1); 
          ysum  = sum(sum(YF/Oth>0.5,1),3);             
          zsum  = shiftdim(sum(sum(YF/Oth>0.5,1),2),1); 
          COM   = [find(cumsum(xsum)/sum(xsum)>0.5,1,'first') ...
                   find(cumsum(ysum)/sum(ysum)>0.5,1,'first') ...
                   find(cumsum(zsum)/sum(zsum)>0.5,1,'first')];
          COMmm = (eye(4) * VF.mat) * [COM';1];
          AffineCOM = eye(4); AffineCOM(13:15) = -COMmm(1:3); 
          %ACvox  = round(inv(Affine * VF.mat) * [ 0; 0; 0; 1]); ACvox = ACvox(1:3)';
          COMvox = round(inv(AffineCOM * VF.mat) * [0;0;0;1]);  COMvox = COMvox(1:3)'; %#ok<MINV,NASGU>
          if ~debug, clear YF xsum ysum zsum Oth; end 



          % load template and remove the sull if the image is skull-stripped
          try 
            VG = spm_vol(Pt1);
          catch
            pause(rand(1))
            VG = spm_vol(Pt1);
          end


          % skull-stripping of the template
          if ppe.affreg.skullstripped
            % print a warning for all users because processing of
            % skull-stripped data is not standard!
            if job.extopts.APP==1 || job.extopts.APP==2
              cat_io_cprintf('warn',[...
                'WARNING: Detected skull-stripping or strongly masked image. Skip APP. \n' ...
                '         Use skull-stripped initial affine registration template and  \n' ...
                '         TPM without head tissue (class 4 and 5)! \n']);
            else
              cat_io_cprintf('warn',[...
                'WARNING: Detected skull-stripping or strongly masked image. \n' ...
                '         Use skull-stripped initial affine registration template and \n' ...
                '         TPM without head tissue (class 4 and 5)! \n']);
            end
            if job.extopts.verb>1
              cat_io_cprintf('warn',sprintf(...
               ['           %0.2f%%%% zeros, %d object(s), %d background region(s) \n' ...
                '           %4.0f cm%s, normalized SD of all tissues %0.2f \n'],...
                ppe.affreg.skullstrippedpara(1:4),char(179),ppe.affreg.skullstrippedpara(5))); 
            end

            % skull-stripping of the template
            VB = spm_vol(Pb);
            [VB2,YB] = cat_vol_imcalc([VG,VB],Pbt,'i1 .* i2',struct('interp',3,'verb',0)); 
            VB2.dat(:,:,:) = eval(sprintf('%s(YB/max(YB(:))*255);',spm_type(VB2.dt))); 
            VB2.pinfo      = repmat([1;0],1,size(YB,3));
            VG             = cat_spm_smoothto8bit(VB2,0.5);
            %VG.dat         = cat_vol_ctype(VG.dat);
            %VG.dt          = [spm_type('UINT8') spm_platform('bigend')];
            clear VB2 YB; 
          end


          % Rescale images so that globals are better conditioned
          VF.pinfo(1:2,:) = VF.pinfo(1:2,:)/spm_global(VF);
          VG.pinfo(1:2,:) = VG.pinfo(1:2,:)/spm_global(VG);


          % APP step 1 rough bias correction 
          % --------------------------------------------------------------
          % Already for the rought initial affine registration a simple  
          % bias corrected and intensity scaled image is required, because
          % high head intensities can disturb the whole process.
          % --------------------------------------------------------------
          %   ds('l2','',vx_vol,Ym, Yt + 2*Ybg,Ysrc/WMth,Ym,60)
          if  job.extopts.APP>2 %app.bias ||  app.msk 
            stime = cat_io_cmd('APP: Rough bias correction'); 
            Ysrc = single(spm_read_vols(spm_vol(obj.image(1).fname)));
            Ysrc(isnan(Ysrc) | isinf(Ysrc)) = min(Ysrc(:));

            [Ym,Yt,Ybg,WMth,bias,Tth,ppe.APPi] = cat_run_job_APP_init(...
              Ysrc,vx_vol,struct('verb',job.extopts.verb,'APPstr',job.opts.biasstr)); 
            bth = min( [ mean(single(Ysrc( Ybg(:)))) - 2*std(single(Ysrc( Ybg(:)))) , ...
                         mean(single(Ysrc(~Ybg(:)))) - 4*std(single(Ysrc(~Ybg(:)))) , ...
                         min(single(Ysrc(~Ybg(:)))) ] );
            if ~debug, clear Ysrc; end 

            % zero background is required for images with high intensity background
            zeroBG = cat_stat_nanmean(Ym(Ybg))<1/3;
            Ymc = Ym; if ~zeroBG, Ymc(Ybg) = 0; bth = 0; end; 
            if ~debug, clear Ym; end 

            % write data to VF
            VF.dt    = [spm_type('UINT8') spm_platform('bigend')];
            VF.pinfo = repmat([1;0],1,size(Ymc,3));
            VF.dat   = cat_vol_ctype(Ymc * max(1/2,Tth.Tmax) * 255); 

            if strcmp('human',job.extopts.species) || 1
              stime = cat_io_cmd('Coarse affine registration:','','',1,stime); 
            end
          else
            if strcmp('human',job.extopts.species) || 1 
              stime = cat_io_cmd('Coarse affine registration'); 
            end
          end

          % smoothing
          VF1   = cat_spm_smoothto8bit(VF,8.0);
          VG1   = cat_spm_smoothto8bit(VG,0.5); 

          % prepare affine parameter 
          aflags     = struct('sep',obj.samp,'regtype','subj','WG',[],'WF',[],'globnorm',1); 
          aflags.sep = max(aflags.sep,max(sqrt(sum(VG(1).mat(1:3,1:3).^2))));
          aflags.sep = max(aflags.sep,max(sqrt(sum(VF(1).mat(1:3,1:3).^2))));



          %% affine registration
          if strcmp('human',job.extopts.species)  || 1
            spm_plot_convergence('Init','Coarse affine registration','Mean squared difference','Iteration');
            warning off 
            try 
              [Affine0, affscale] = cat_spm_affreg(VG1, VF1, aflags,Affine, 1); Affine = Affine0; 
            catch
              Affine0 = eye(4); affscale = 1; cat_io_cmd(' ','','',1);
            end
            warning on


            %  If this registration failed and the AC is somewhere outside 
            %  than use an initial affine registion that use the COM as AC. 
            %  Use the more successfull mapping
            mat0 = spm_imatrix(Affine0);
            ACvox  = round(inv(Affine * VF.mat) * [ 0; 0; 0; 1]); ACvox = ACvox(1:3)'; %#ok<MINV>
            if 0 && ...
              any(any(isnan(Affine0(1:3,:)))) || affscale<0.5 || affscale>3 || (all(all(roundx(Affine0,6) == eye(4))) &&  affscale == 1) || ...
              (max(mat0(7:9))/min(mat0(7:9)))>1.2 || any(ACvox<VF.dim/6) || any(ACvox>VF.dim - VF.dim/6)

              % do affreg
              warning off 
              try 
                [AffineCOM2, affscaleCOM]  = cat_spm_affreg(VG1, VF1, aflags, AffineCOM, 1);
              catch
                AffineCOM2 = eye(4); affscaleCOM = 1; 
              end
              warning on

              matCOM  = spm_imatrix(AffineCOM2);
              COMvox2 = round(inv(AffineCOM2 * VF.mat) * [0;0;0;1]);  COMvox2 = COMvox2(1:3)'; %#ok<MINV>

              % check affreg and use the affine registration with less different scaling factors 
              if (( (max(matCOM(7:9))/min(matCOM(7:9))) < (max(mat0(7:9))/min(mat0(7:9)))*1.05 ) && ...
                 (  abs(affscaleCOM-1) < abs(affscale-1)*1.05 )) || ...
                 (any(COMvox2>VF.dim/6) && any(COMvox2<VF.dim - VF.dim/6))
                Affine = AffineCOM2; affscale = affscaleCOM; 

                % give another (silent) warning
                if job.extopts.verb>1
                  cat_io_cprintf('warn','\nInitial registration failed use center of mass as AC! \n');
                  cat_io_cmd(' ','','',1);  
                end
              end
            end

            % final check of the new affine matrix 
            if any(any(isnan(Affine(1:3,:)))) || affscale<0.5 || affscale>3 || (all(all(roundx(Affine,6) == eye(4))) &&  affscale == 1)
              stime  = cat_io_cmd('Coarse affine registration failed. Try fine affine registration.','','',1,stime);
              cat_io_cmd(' ','','',1);  Affine = eye(4); affscale = 1; cat_io_cmd(' ','','',1);
            end

          else
            % no affine registration 
            Affine = eye(4); affscale = 1; Affine0 = eye(4); %#ok<NASGU> % Affine 0 for debuging
          end




          %% APP step 2 - brainmasking and second tissue separated bias correction  
          %  ---------------------------------------------------------
          %  The second part of APP maps a brainmask to native space and 
          %  refines it by morphologic operations and region-growing to
          %  adapt for worse initial affine alignments. It is important
          %  that the mask covers the whole brain, whereas additional
          %  masked head is here less problematic.
          %  ---------------------------------------------------------
          %    ds('l2','',vx_vol,Ym,Yb,Ym,Yp0,90)
          aflags.sep = obj.samp/2; 
          aflags.sep = max(aflags.sep,max(sqrt(sum(VG(1).mat(1:3,1:3).^2))));
          aflags.sep = max(aflags.sep,max(sqrt(sum(VF(1).mat(1:3,1:3).^2))));

          %% apply (first affine) registration on the default brain mask
          VFa = VF; if strcmp(job.extopts.species,'human') || 1, VFa.mat = Affine * VF.mat; else Affine = eye(4); affscale = 1; end
          if isfield(VFa,'dat'), VFa = rmfield(VFa,'dat'); end
          [Vmsk,Yb] = cat_vol_imcalc([VFa,spm_vol(Pb)],Pbt,'i2',struct('interp',3,'verb',0)); 
          if exist('Ybg','var'), Yb = Yb>0.5 & ~Ybg; end

          if exist('Ymc','var') && (sum(Yb(:))==0 || mean(Ymc(Yb(:)>0))<0.5 || mean(Ymc(Yb(:)>1.05))<0.5) 
            %%
            [Ymr,Ytr,resT3] = cat_vol_resize({Ymc,single(Yt)},'reduceV',vx_vol,min(1.2,cat_stat_nanmean(vx_vol)*2),64,'meanm'); 
            Ybr  = cat_vol_morph(Ytr>0.5 & Ymr>0.7 & Ymr<1.2,'o',1); 
            Ybsr = cat_vol_smooth3X(Ybr,8);  
            Ybsr = Ybsr./max(Ybsr(:)); 
            Ybr  = cat_vol_morph(Ytr>0.3 & Ymr>0.7 & Ymr<1.2 & Ybsr>0.1,'o',1); 
            Ybr  = cat_vol_morph(Ybr,'lc',6); 
            Yb   = cat_vol_resize(Ybr,'dereduceV',resT3); 
            Affine = eye(4); affscale = 1;
          end

          if job.extopts.APP>3
            stime = cat_io_cmd('APP: Fine bias correction and skull-stripping:','','',1,stime); 

            % fine APP
            [Ym,Yp0,Yb] = cat_run_job_APP_final(single(obj.image.private.dat(:,:,:)+0),...
                Ymc,Yb,Ybg,vx_vol,job.extopts.gcutstr,job.extopts.verb);
            stime = cat_io_cmd('Affine registration:','','',1,stime); 

            % zero background is required for images with high intensity background
            Ymc = Ym; if ~zeroBG, Ymc(Ybg) = 0; end
            if ~debug, clear Ym; end

            VF.dat(:,:,:) = cat_vol_ctype(Ymc * max(1/2,Tth.Tmax) * 255); 

          end
          if strcmp('human',job.extopts.species) || 1
            stime = cat_io_cmd('Affine registration','','',1,stime); 
            spm_plot_convergence('Init','Affine registration','Mean squared difference','Iteration');

            % smooth data
            VF1 = cat_spm_smoothto8bit(VF,aflags.sep);
            VG1 = cat_spm_smoothto8bit(VG,0.5); 

            % fine affine registration 
            warning off
            try
              [Affine1,affscale1] = spm_affreg(VG1, VF1, aflags, Affine, affscale); 
            catch
              Affine1 = eye(4,4); affscale1 = 1;
            end
            warning on

            % check results and reset affine matrix 
            if any(any(isnan(Affine1(1:3,:)))) || affscale1<0.5 || affscale1>3,
              Affine1 = eye(4,4); affscale1 = 1;
            end
            Affine = Affine1; 

            if ~debug, clear VG VG1 VF1; end

            ppe.affreg.Affine1 = Affine; ppe.affreg.mat1 = spm_imatrix(Affine); ppe.affreg.affscale1 = affscale1;
          end


          %% APP for spm_maff8
          %  optimize intensity range
          %  we have to rewrite the image, because SPM reads it again 
          if job.extopts.APP>3 %&& (app.aff>0 || app.bias>0 || app.msk>0)
              % add and write bias corrected (, skull-stripped) image
              Ymc2 = single(max(bth,min( max( WMth + max(3*diff([bth,WMth])) , ...
                WMth / max(1/2 - 3/8*Tth.inverse,Tth.Tmax) ) ,Ymc * WMth))); % use the bias corrected image
              if ~debug, clear Ymc; end

              % hard masking
              %if ~zeroBG, Ymc2 = Ymc2 .* (~Ybg); end % 20161229 - simple to do this than to change all things in cat_main
%              if job.extopts.APP==5, Ymc2 = Ymc2 .* cat_vol_morph(cat_vol_morph(Yb,'d',1),'lc',1); end

              % set variable and write image
              %obj.image.dat   = Ymc2;  
              %obj.image.dt    = [spm_type('FLOAT32') spm_platform('bigend')];
              %obj.image.pinfo = repmat([1;0],1,size(Ymc2,3));
              spm_write_vol(obj.image,Ymc2); 
          end
        

        
          %%
          if ppe.affreg.skullstripped  
            %% update number of SPM gaussian classes 
            Ybg = 1 - spm_read_vols(obj.tpm.V(1)) - spm_read_vols(obj.tpm.V(2)) - spm_read_vols(obj.tpm.V(3));
            obj.tpm.V(4).dat = Ybg; 
            obj.tpm.dat{4}   = Ybg; 
            obj.tpm.V(4).pinfo = repmat([1;0],1,size(Ybg,3));
            obj.tpm.dat(5:6) = []; 
            obj.tpm.V(5:6)   = []; 
            obj.tpm.bg1(4)   = obj.tpm.bg1(6);
            obj.tpm.bg2(4)   = obj.tpm.bg1(6);
            obj.tpm.bg1(5:6) = [];
            obj.tpm.bg2(5:6) = [];

            %job.opts.ngaus = [job.opts.ngaus(1:3) 1];
            job.opts.ngaus = 3*ones(4,1); % this is more save
            obj.lkp        = [];
            for k=1:numel(job.opts.ngaus)
              job.tissue(k).ngaus = job.opts.ngaus(k);
              obj.lkp = [obj.lkp ones(1,job.tissue(k).ngaus)*k];
            end
          end
          %% guarantee SPM image requirements of gaussian distributed data
          %  Add noise to regions with low intensity changes that are mosty 
          %  caused by skull-stripping, defacing, or simulated data and can
          %  lead to problems in SPM tissue peak estimation. 
          if 1
            %%

            Ysrc  = single(spm_read_vols(obj.image));

            Yg    = cat_vol_grad(Ysrc,vx_vol);
            Ygnth = cat_stat_nanmean( min(0.3,Yg(:)./max(eps,Ysrc(:))) );
            gno   = Ysrc( ( Yg(:)./max(eps,Ysrc(:)) )<Ygnth/2 & Ysrc(:)>cat_stat_nanmean(Ysrc(:)) ); 
            Tthn  = mean( Ysrc( ( Yg(:)./max(eps,Ysrc(:)) )<Ygnth/2 & ...
                      Ysrc(:)>(cat_stat_nanmean(gno) - 1*cat_stat_nanstd(gno)) & ...
                      Ysrc(:)<(cat_stat_nanmean(gno) + 4*cat_stat_nanstd(gno))) ); 
            Ysrcn = max( min(Ysrc(:)) , Ysrc + ( randn(size(Ysrc))*Tthn*0.005 + (rand(size(Ysrc))-0.5)*Tthn*0.01 ) .*  ...
              cat_vol_smooth3X( Ysrc/Tthn<0.01 | Yg/Tthn<0.01 | Yg./max(eps,Ysrc)<0.01 , 4 ));

            [h,i] = hist(Ysrcn(Ysrc~=0 & ~isnan(Ysrc)),1000); maxth = cumsum(h);
            idl = find(maxth/max(maxth)>0.02,1,'first'); 
            idh = find(maxth/max(maxth)>0.98,1,'first'); 
            iil = i(idl) - diff([idl idh])*0.5; 
            iih = i(idh) + diff([idl idh])*0.5; 
            Ysrcn = max( iil , min( iih , Ysrcn )); 

            % should be used above .. however SPM don't like it ... 
            if ~strcmp(job.extopts.species,'human') && exist('Yb','var')
              %Ysrcn = Ysrcn .* cat_vol_morph(Yb,'d',3); 
              %ppe.affreg.skullstripped = 1;
            end


            %%
            spm_write_vol(obj.image,Ysrcn); 
            if ~debug, clear Ysrc Yg Ysrcn; end
          end


        
          %%  Fine Affine Registration with 3 mm sampling distance
          %  This does not work for non human (or very small brains)
          stime = cat_io_cmd('SPM preprocessing 1 (estimate):','','',1,stime); 
          if strcmp('human',job.extopts.species) || 1
            % sampling >6 mm seams to fail in some cases, but higher sampling (<3 mm) did not improve the result 
            % also the variation of smoothing did not work well in all cases and lower smoothing (<=4 mm) lead to problems in some cases 
            fwhm = [0 24 12 8]; samp = [0 6 6 6]; % entry 1 is the initial affine registration 
            ppe.affreg.Affine2   = cell(size(fwhm));  ppe.affreg.Affine2{1}   = Affine1; 
            ppe.affreg.affscale2 = cell(size(fwhm));  ppe.affreg.affscale2{1} = ppe.affreg.affscale1; 
            ppe.affreg.mat2      = cell(size(fwhm));  ppe.affreg.mat2{1}      = ppe.affreg.mat1; 
            ppe.affreg.useaff2   = zeros(size(fwhm)); ppe.affreg.useaff2(1)   = 1; 
            ppe.affreg.stopiter  = zeros(size(fwhm));
            ppe.affreg.affll     = inf(size(fwhm));
            ppe.affreg.sc1{1}    = (max(ppe.affreg.mat2{1}(7:9)) / min(ppe.affreg.mat2{1}(7:9)));
            ppe.affreg.sc2{1}    = abs(ppe.affreg.affscale2{1} - 1);
            %
            for fwhmi=2:numel(fwhm)
              %%
              spm_plot_convergence('Init','Fine affine registration','Mean squared difference','Iteration');
              regtime = clock;

              % spm_maff8 registration
              warning off
              [ppe.affreg.Affine2{fwhmi}, ppe.affreg.affll(fwhmi), ppe.affreg.affh{fwhmi}] = ...
                spm_maff8( obj.image(1) , samp(fwhmi) , fwhm(fwhmi) , obj.tpm , ppe.affreg.Affine2{fwhmi} , job.opts.affreg );
              warning on  

              ppe.affreg.mat2{fwhmi}       = spm_imatrix(ppe.affreg.Affine2{fwhmi}); 
              ppe.affreg.affscale2{fwhmi}  = mean(ppe.affreg.mat2{fwhmi}(7:9)); 

              % check registration and reset values to previous results
              ppe.affreg.sc1{fwhmi} = (max(ppe.affreg.mat2{fwhmi}(7:9)) / min(ppe.affreg.mat2{fwhmi}(7:9)));
              ppe.affreg.sc2{fwhmi} = abs(ppe.affreg.affscale2{fwhmi} - 1);
              ppe.affreg.useaff2(fwhmi) = fwhmi<4 || ...
                (ppe.affreg.sc1{fwhmi} < ppe.affreg.sc1{fwhmi-1}*1.05 || ...            % check scaling properties
                 ppe.affreg.sc2{fwhmi} < ppe.affreg.sc2{fwhmi-1}*1.05) && ...            % ckeck total scaling
                ppe.affreg.affll(fwhmi) < ppe.affreg.affll(fwhmi-1)*1.05 && ...
                ppe.affreg.affll(fwhmi) > 0.5 && ppe.affreg.affll(fwhmi) < 1.5 && ...  % check SPM convergence criteria
                ppe.affreg.affscale2{fwhmi}>0.5 && ppe.affreg.affscale2{fwhmi}<3;      % check principle scaling range
              ppe.affreg.stopiter(fwhmi) = fwhmi>3 && (...
                ppe.affreg.useaff2(fwhmi)==0 || ... 
                ppe.affreg.sc1{fwhmi}/ppe.affreg.sc1{fwhmi-1}>1.05 || ... % stop if values get worse
                ppe.affreg.sc1{fwhmi}/ppe.affreg.sc1{fwhmi-1}>1.05 || ...
                ppe.affreg.affll(fwhmi)/ppe.affreg.affll(fwhmi-1)>1.05 || ... % stop for low changes
                abs(ppe.affreg.sc1{fwhmi} - ppe.affreg.sc1{fwhmi-1})<0.01 || ...
                abs(ppe.affreg.sc1{fwhmi} - ppe.affreg.sc1{fwhmi-1})<0.01 || ...
                ((ppe.affreg.affll(fwhmi) - ppe.affreg.affll(fwhmi-1))<0.01 && fwhmi>2));

              %% some information in case of debugging
              if 0 || debug || job.extopts.verb > 2
                if fwhmi==2, fprintf('\n'); end; 
                fprintf('  sc: %5.3f >> %5.3f; sc2: %5.3f >> %5.3f; conv: %5.3f > %5.3f, time: %3.0fs - use: %d, stop: %d\n',...
                  (max(ppe.affreg.mat2{fwhmi-1}(7:9)) / min(ppe.affreg.mat2{fwhmi-1}(7:9))), ...
                  (max(ppe.affreg.mat2{fwhmi}(7:9))   / min(ppe.affreg.mat2{fwhmi}(7:9))), ...
                  abs(ppe.affreg.affscale2{fwhmi-1} - 1), abs(ppe.affreg.affscale2{fwhmi} - 1), ...
                  ppe.affreg.affll(fwhmi-1), ppe.affreg.affll(fwhmi) , etime(clock,regtime) , ...
                  ppe.affreg.useaff2(fwhmi) ,ppe.affreg.stopiter(fwhmi) ); 
                if ppe.affreg.stopiter(fwhmi)
                  cat_io_cmd(' ','','',1); 
                end
              end

              if ppe.affreg.useaff2(fwhmi)==0
                ppe.affreg.Affine2{fwhmi}   = ppe.affreg.Affine2{fwhmi-1};
                ppe.affreg.mat2{fwhmi}      = ppe.affreg.mat2{fwhmi-1}; 
                ppe.affreg.affll(fwhmi)     = ppe.affreg.affll(fwhmi-1); 
                ppe.affreg.affscale2{fwhmi} = ppe.affreg.affscale2{fwhmi-1}; 
              end
               Affine = ppe.affreg.Affine2{fwhmi}; 

              %% stop iteration if the results bring no advante 
              if ppe.affreg.stopiter(fwhmi), break; end 
            end
          end
          obj.Affine = Affine;


          if 0
            %% just for debugging!
            if ~exist('Ym','var'), Vm = cat_spm_smoothto8bit(VF,0.5); Ysrc = Vm.dat; WMth=255; Ym = single(Ysrc)/WMth*2; end
            VFa = VF; if isfield(VFa,'dat'), VFa = rmfield(VFa,'dat'); end
            if 1 
              if 0 & isfield(ppe.affreg,'Affine2') %&& ~isempty(ppe.affreg.Affine2{end})
                VFa.mat = ppe.affreg.Affine2{2} * VF.mat; 
              else
                VFa.mat = Affine * VF.mat;
              end
            else % old registration
              VFa.mat = Affine0 * VF.mat;
            end
            [Vmsk,Yb] = cat_vol_imcalc([VFa,spm_vol(Pb)],Pbt,'i2',struct('interp',3,'verb',0)); 
            %[Ysrcr,Ybr] = cat_vol_resize({Ysrc,Yb},'reduceV',vx_vol,2,32,'meanm'); Tth0 = kmeans3D(Ysrcr(Ybr(:)>0.5),3);
            %Ym = Ysrc/Tth0(3); %cat_stat_nanmean(Ysrc(Ysrc(:)>cat_stat_nanmean(Ysrc(:)))); 
            if ~exist('WMth','var'), WMth = cat_stat_nanmean(Ysrc(Yb(:)>0.5 & Ysrc(:)>mean(Ysrc(Yb(:)>0.5)))); end
            AC = round(inv(Affine * VF.mat) * [ 0; 0; 0; 1]);
            if 1
              %%
              switch 2
                case 1, ds('l2','',vx_vol,Ym,Yb>0.5,Ysrc/WMth,Ym,AC(3))
                case 2, ds('l2','a',vx_vol,Ym,Yb>0.5,Ysrc/WMth,Ym,AC(2))
                case 3, ds('l2','m',vx_vol,Ym,Yb>0.5,Ysrc/WMth,Ym,AC(1))
              end
            end
          end
          cat_err_res.obj = obj; 
        else
          stime = cat_io_cmd('SPM preprocessing 1 (estimate):','','',1,stime); 
          obj.Affine = Affine; 
        end
        
        
        % set original non-bias corrected image in case of APP_init
        if job.extopts.APP==3
          obj.image = spm_vol(obj.image);
        end

if job.extopts.APP==2 && ~strcmp(job.extopts.species,'human')        
  obj.biasreg  = 0.001;
  obj.biasfwhm = 50;
end
        


        
        
        
        %% SPM preprocessing 1
        %  ds('l2','a',0.5,Ym,Ybg,Ym,Ym,140);
        %  ds('l2','a',0.5,Ysrc/WMth,Yb,Ysrc/WMth,Yb,140);
        warning off 
        try 
          res = spm_preproc8(obj);
        catch
          % only in default mode
          if job.extopts.expertgui==0 
            job2          = job;
            job2.channel(1).vols{subj} = job.channel(1).vols0{subj};
            if job.extopts.APP<1
              cat_io_cprintf('err','\n  Failed, try APP=1! \n'); 
              job2.extopts.APP = 1; 
            elseif job.extopts.APP<4
              cat_io_cprintf('err','\n  Failed, try APP=4! \n'); 
              job2.extopts.APP = 4; 
            end
            cat_run_job(job2,tpm,subj); 
            return
          end
          if ppe.affreg.skullstripped  
            obj.tpm = tpm; 
          end
          
          
          if job.extopts.NCstr || any( (vx_vol ~= vx_voli) ) || ~strcmp(job.extopts.species,'human')
            [pp,ff,ee] = spm_fileparts(job.channel(1).vols{subj});
            delete(fullfile(pp,[ff,ee]));
            error('CAT:cat_run_job:spm_preproc8','Error in spm_preproc8. Check image and orientation. \n');
          end
        end
        
        warning on 
        if exist('ppe','var'), res.ppe = ppe; end

        if job.extopts.experimental
            % save information for debuging and OS test
            [pth,nam] = spm_fileparts(job.channel(1).vols0{subj}); 
            tmpmat = fullfile(pth,reportfolder,sprintf('%s_%s_%s.mat',nam,'runjob','postpreproc8')); 
            save(tmpmat,'obj','res','Affine','Affine0','Affine1','Affine2');     
        end 
        cat_err_res.res = res;   

        fprintf('%5.0fs\n',etime(clock,stime));   


        %% check contrast  
        clsint = @(x) round( sum(res.mn(res.lkp==x) .* res.mg(res.lkp==x)') * 10^5)/10^5;
        Tgw = [cat_stat_nanmean(res.mn(res.lkp==1)) cat_stat_nanmean(res.mn(res.lkp==2))]; 
        Tth = [
          ... min(res.mn(res.lkp==6 & res.mg'>0.3)) ... % bg; ignore the background, because of MP2RGAGE, R1, and MT weighted images  
          max( min( clsint(3) ,  max(Tgw)+abs(diff(Tgw))) , min(Tgw)-abs(diff(Tgw)) ) ... % csf with limit for T2!
          clsint(1) ... gm
          clsint(2) ... wm 
        ];
        
        % save data for error report
        %if isfield(obj,'msk'), res.msk = obj.msk; end
        res.Tth = Tth; 
        cat_err_res.res = res;   
        
        % inactive preprocessing of inverse images (PD/T2) 
        if job.extopts.INV==0 && any(diff(Tth)<=0)
          error('CAT:cat_main:BadImageProperties', ...
          ['CAT12 is designed to work only on highres T1 images.\n' ...
           'T2/PD preprocessing can be forced on your own risk by setting \n' ...
           '"cat12.extopts.INV=1" in the cat default file. If this was a highres \n' ...
           'T1 image then the initial segmentation might be failed, probably \n' ...
           'because of alignment problems (please check image orientation).']);    
        end
        
        
    end
    
    
    %% call main processing
    res.stime  = stime;
    res.catlog = catlog; 
    res.image0 = spm_vol(job.channel(1).vols0{subj}); 
    cat_main(res,obj.tpm,job);
    
    % delete denoised/interpolated image
    [pp,ff,ee] = spm_fileparts(job.channel(1).vols{subj});
    if exist(fullfile(pp,[ff,ee]),'file'); 
      delete(fullfile(pp,[ff,ee]));
    end
%%
return
%=======================================================================
function r = roundx(r,rf)
  r(:) = round(r(:) * rf) / rf;
return
%=======================================================================
