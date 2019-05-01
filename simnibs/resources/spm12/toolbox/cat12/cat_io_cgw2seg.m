function varargout=cat_io_cgw2seg(c,g,w,opt)
% ______________________________________________________________________
% convert SPM or FSL probability maps to PVE label files
%
%   varargout=cat_io_cgw2seg(c,g,w,mode,d)
%
%   mode .. {'SPM','FLS'}
%   d    .. delete old files [0|1], default = 0;
%
% ______________________________________________________________________
% $Revision: 1139 $  $Date: 2017-06-12 21:54:51 +0200 (Mo, 12 Jun 2017) $

  if ~exist('opt','var'); opt = struct(); end
  
  def.verb   = 1;
  def.mode   = 'SPM';
  def.delete = 0;
  def.lazy   = 1; 
  def.fsd    = '/Volumes/vbmDB/MRData/vbm12tst/results/deffiles/fs6';
  def.fsh    = '/Applications/freesurfer'; 
  def.fss    = '/Applications/freesurfer/subjects'; 
  opt        = cat_io_checkinopt(opt,def);
  opt.fsc    = ['export FREESURFER_HOME=' def.fsh '; source $FREESURFER_HOME/SetUpFreeSurfer.sh; ']; 
  
  if opt.verb, spm('fnbanner','cat_io_cgw2seg'); end
  % if there is no input, select files
  if ~exist('c','var') || ~exist('g','var') || ~exist('w','var') || isempty(c) || isempty(g) || isempty(w)
    if ~isfield(opt,'mode') || isempty(opt.mode)
      c = spm_select(Inf  ,'image','Select CSF files'); n = size(c,1); [wd,ff]=spm_fileparts(c(1,:));
      if strcmp(ff(1:2),'c1') || strcmp(ff(1:2),'c2') || strcmp(ff(1:2),'c3') 
        g=c; w=c; 
        for i=1:n
          g(i,numel(wd) + 3)='1'; 
          w(i,numel(wd) + 3)='2'; 
          c(i,numel(wd) + 3)='3'; 
        end
      else
        g = spm_select([n n],'image','Select GM  files',{},wd);
        w = spm_select([n n],'image','Select WM  files',{},wd);
      end
    else
      switch opt.mode
        case 'FS'
          % select subject directories
          fsdirs = cellstr(spm_select(Inf,'dir','Select subject dirs',{},opt.fss));
          % get result directory
          if ~exist(opt.fsd,'dir')
            ndir  = cellstr(spm_select(1  ,'dir','Select output directory')); 
          else
            ndir  = opt.fsd;
          end
          % find segmentation files
          pamgz = {}; p = {};
          for i=1:numel(fsdirs)
            [tmp,name] = spm_fileparts(fsdirs{i}); 
            file = fullfile(fsdirs{i},'mri','aseg.auto_noCCseg.mgz'); 
            if exist(file,'file')
              % set output filename
              [h,f]  = fileparts(file);
              [pp,ff] = spm_fileparts(spm_fileparts(h)); 
              f = cat_vol_findfiles(opt.fsd,[ff '.nii']); 
              if ~isempty(f)
                cat_io_cprintf([0 0.5 0],'Found segmentation of "%s". \n',name)
                pamgz{end+1} = file; 
                p{end+1} = fullfile(fileparts(f{1}),['p0' ff '.nii']);
              else
                cat_io_cprintf([0.5 0.5 0],'Found segmentation of "%s", but not the original input nifti. \n',name)
                pamgz{end+1} = file; 
                p{end+1} = fullfile(opt.fsd,['p0' ff '.nii']);
              end
            else
              cat_io_cprintf([0.5 0 0],'Found no segmentation of "%s". \n',name)
            end
          end
          
          %% convert mgz to nifti
          for i=1:numel(pamgz) 
            pa{i} = [pamgz{i}(1:end-3) 'nii']; c{i} = pa{i}; g{i} = pa{i}; w{i} = pa{i};
            psmgz{i} = strrep(pamgz{i},'aseg.auto_noCCseg.mgz','brainmask.auto.mgz'); 
            ps{i} = [psmgz{i}(1:end-3) 'nii'];
            pomgz{i} = strrep(pamgz{i},'aseg.auto_noCCseg.mgz',fullfile('orig','001.mgz')); 
            po{i} = [pomgz{i}(1:end-3) 'nii'];
            if ~opt.lazy || ~exist(pa{i},'file') || ~exist(ps{i},'file') 
              [tmp,name] = spm_fileparts(p{i}); fprintf('Convert mgz to nii of "%s". \n',name)
              [SR,SR] = system(sprintf('%smri_convert -i %s -o %s -ot nii',opt.fsc,pamgz{i},pa{i}));
              [SR,SR] = system(sprintf('%smri_convert -i %s -o %s -ot nii',opt.fsc,psmgz{i},ps{i}));
              if strcmp(spm_fileparts(p{i}),opt.fsd)  
                [SR,SR] = system(sprintf('%smri_convert -i %s -o %s -ot nii',opt.fsc,pomgz{i},po{i}));
              end
            end
          end
       
         
          
        case 'FSLs'
          c2 = cellstr(spm_select(Inf  ,'any','Select CSF files',{},'','T1_fast_pve_0.nii.*')); 
          if isempty(c2) || isempty(c2{1}), return; end
          for i=1:numel(c2), if strcmp(c2{i}(end-2:end),'nii'), c2{i} = [c2{i} '.gz']; end; end
          g2=c2; w2=c2; for i=1:numel(c2), g2{i}(end-7)='1'; w2{i}(end-7)='2'; end
          c=c2; g=g2; w=w2; for i=1:numel(c), c{i}(end-2:end)=[]; g{i}(end-2:end)=[]; w{i}(end-2:end)=[]; end
          gzfiles = [c2;g2;w2]; 
          for i=1:numel(gzfiles),
            if exist(gzfiles{i},'file'), system(sprintf('gunzip %s',gzfiles{i})); end
          end
        case 'FSL', 
          c = cellstr(spm_select(Inf  ,'image','Select CSF files',{},'','*seg_0.*')); 
          if isempty(c) || isempty(c{1}), return; end
          g=c; w=c; for i=1:numel(c), g{i}(end-6)='1'; w{i}(end-6)='2'; end
        case 'SPM', 
          g = cellstr(spm_select(Inf  ,'image','Select GM files',{},'','^c1.*'));       
          if isempty(g) || isempty(g{1}), return; end
          c=g; w=g; p=g; 
          for i=1:numel(g), 
            [h,f,e,nn]=spm_fileparts(g{i}); e=e(1:4); 
            g{i}=fullfile(h,[f,e,nn]); f(2)='3'; 
            c{i}=fullfile(h,[f,e,nn]); f(2)='2'; 
            w{i}=fullfile(h,[f,e,nn]); f(1:2)='p0'; 
            p{i}=fullfile(h,[f,e]); 
          end
        otherwise 
          c = spm_select(Inf  ,'image','Select CSF files'); n = size(c,1); wd=spm_fileparts(c(1,:));
          g = spm_select([n n],'image','Select GM  files',{},wd);
          w = spm_select([n n],'image','Select WM  files',{},wd);
      end
    end
  end
  
  % if the input is given by char, convert it to cellstr
  if isa(c,'char'), c=cellstr(c); g=cellstr(g); w=cellstr(w); end

  % if the input is give by a c,g,w matrix that only one loop
  if isa(c,'cell'), n=numel(c); else n=1; end

  % remove , from filenames (spm-image number)
  if exist('c','var')
    for ci=1:numel(c),
      [pp,ff,ee] = spm_fileparts(c{ci}); c{ci}=fullfile(pp,[ff ee]);
      [pp,ff,ee] = spm_fileparts(g{ci}); g{ci}=fullfile(pp,[ff ee]);
      [pp,ff,ee] = spm_fileparts(w{ci}); w{ci}=fullfile(pp,[ff ee]);
    end
  end
  
  
  spm_progress_bar('Init',n,'Filtering','Volumes Complete');
  for i=1:n
    if opt.verb, fprintf('%s: ',c{i}); end
    tic
    if isa(c,'cell')
      if opt.lazy && exist('p','var') && exist(p{i},'file') 
        fprintf(' already exist (lazy mode). \n');
        continue
      else
        if exist('pa','var')
          ha = spm_vol(pa{i}); A=uint8(spm_read_vols(ha)); hc=ha;
          hs = spm_vol(ps{i}); S=single(spm_read_vols(hs));
         %%
          C  = single(A== 4 | A== 5 | A==14 | A==15 | A==72 | A==24 | ...
                      A==43 | A==44 | A==64 | A==63 | A==31);
          G  = single(A== 3 | A== 8 | A== 9 | A==10 | A==11 | A==12 | A==13 | A==17 | A==18 | A==30 | A==26 | ...
                      A==42 | A==47 | A==48 | A==49 | A==50 | A==51 | A==52 | A==53 | A==54 | A==57 | A==58 ); 
          W  = single(A==1  | A== 2 | A== 7 | A==16 | A==28 | ...
                      A==40 | A==41 | A==46 | A==55 | A==60 );
          H  = single(A==77 | A==78 | A==79); 
          W  = W+H; 
          C  = C | (A==0 & S~=0);  
          
          %SEG = C + 2*G + 3*W; ds('l2','a',1,S/median(S(S(:)>0)),single(A)/50,SEG/3,single(A)/50,120)
          
        elseif exist(c{i},'file') && exist(g{i},'file') && exist(w{i},'file')
          try
            hc=spm_vol(c{i}); C=spm_read_vols(hc); 
            hg=spm_vol(g{i}); G=spm_read_vols(hg); 
            hw=spm_vol(w{i}); W=spm_read_vols(hw); 
          catch
            fprintf(1,'ERROR:cat_io_cgw2seg - incorrect file(s) [CSF=%d;GM=%d;WM=%d]\n',...
            ~exist(c{i},'file'),~exist(g{i},'file'),~exist(w{i},'file'));
            continue
          end
          %[h,f]  = fileparts(hc.fname);
          %if      strfind({'c1','c2','c3'},f),          mode='SPM'; 
          %elseif  strfind({'_seg0','_seg1','_seg2'},f), mode='FSL';
          %else                                          mode='';
          %end
        else 
          fprintf(1,'ERROR:cat_io_cgw2seg - miss file(s) [CSF=%d;GM=%d;WM=%d]\n',...
            ~exist(c{i},'file'),~exist(g{i},'file'),~exist(w{i},'file'));
          continue
        end
      end
    end

    SEG = C + 2*G + 3*W;
    if max(C(:))>1 || max(G(:))>1 || max(W(:))>1, SEG = SEG/max(W(:)); end

    if exist('hc','var')
      %%
      [h,f]  = fileparts(hc.fname);
      if opt.delete==1 && ~(strcmp(opt.mode,'FS') && ~isempty(ndir)) 
        switch hc.fname(end-2:end)
          case 'nii'
            delete(hc.fname); delete(hg.fname); delete(hw.fname);
          case 'img'
            delete([hc.fname(1:end-3) 'hdr']); delete([hc.fname(1:end-3) 'img']);
            delete([hg.fname(1:end-3) 'hdr']); delete([hg.fname(1:end-3) 'img']);
            delete([hw.fname(1:end-3) 'hdr']); delete([hw.fname(1:end-3) 'img']);
        end
      end
      switch opt.mode
        case 'FS',   hc.fname = fullfile(h,['p0' f(1:end)   '.nii']);
        case 'SPM',  hc.fname = fullfile(h,['p0' f(3:end)   '.nii']);
        case 'FSL',  hc.fname = fullfile(h,['p0' f(1:end-6) '.nii']);
        case 'FSLs', 
          [h2,f2] = spm_fileparts(h); 
          hc.fname = fullfile(h2,['p0' f2 '.nii']);
        otherwise,   hc.fname = fullfile(h,['p0' f          '.nii']);
      end
      %hc.dt(1) = 16;
      hc.dt(1) = spm_type('uint8');
      hc.pinfo = [3/255;0;1];
      hc.descript = 'label map';
      if exist(hc.fname,'file'),delete(hc.fname); end
      spm_write_vol(hc,SEG);

      %%
      if strcmp(opt.mode,'FS') && ~isempty(ndir) && ~(opt.lazy && exist('p','var') && exist(p{i},'file')) 
        %% find native file 
        [pp,ff] = spm_fileparts(spm_fileparts(h)); 
        ofile   = cat_vol_findfiles(ndir,[ff '.nii']);
        if isempty(ofile), cat_io_cprintf([0.5 0 0], 'No original file!\n'); continue; end 
          
        clear Vi Vo; 
        
        Vi = spm_vol(char([ofile;{hc.fname}])); 
        Vo = Vi(1); Vo.fname = p{i}; Vo.pinfo = [3/255;0;1]; Vo.dt = hc.dt; Vo.descrip = 'label map';
        
        % map to origal space
        cat_vol_imcalc(Vi,Vo,'single(i2)',struct('interp',1)); 
        cat_io_cprintf([0 0.5 0],'done');
      end
      
      if nargout>0, varargout{1}{n} = hc.fname; end
    end
    spm_progress_bar('Set',i); if opt.verb, fprintf('%4.0f\n',toc); end
  end
  spm_progress_bar('Clear');
  if opt.verb, fprintf('%-40s: %30s\n','Completed',spm('time')); end
end
