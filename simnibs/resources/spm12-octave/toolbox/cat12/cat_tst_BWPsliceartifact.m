%% create BWP slice artifact phantom
%  ---------------------------------------------------------------------
%  Use image without bias for better testing, as far as this correction
%  will also reduce the main inhomogeneity and not only the slice artifact.

%BWP image
Pp0  = '/Volumes/MyBook/MRData/cat_tst/+RAW/BWP_Collins/SEG/p0BWPC_HC.nii';
P    = {
%  '/Volumes/4TBWD/MRData/database/BWP/results/BWPC_NIR/T1/BWPC_HC_T1_pn1_rf00000_vx100x100x100.nii'
  '/Volumes/4TBWD/MRData/database/BWP/results/BWPC_NIR/T1/BWPC_HC_T1_pn3_rf040pA_vx100x100x100.nii'
  '/Volumes/4TBWD/MRData/database/BWP/results/BWPC_NIR/T1/BWPC_HC_T1_pn7_rf080pA_vx100x100x100.nii'
%  '/Volumes/4TBWD/MRData/database/BWP/results/BWPC_NIR/T1/BWPC_HC_T1_pn3_rf080pA_vx100x100x100.nii'
%  '/Volumes/4TBWD/MRData/database/BWP/results/BWPC_NIR/T1/BWPC_HC_T1_pn7_rf040pA_vx100x100x100.nii'
  };
rdir = '/Volumes/4TBWD/MRData/database/BWP/sliceart';
if ~exist(rdir,'dir'), mkdir(rdir); end

% parametersets
art.intraSlAmp  = 0.30:0.30:0.90; 
art.interSlAmp  = 0.20:0.20:0.60; 
art.intraSlFreq = [2 4 8];
art.interSlFreq = [0 2 8]; 

for pi=1:numel(P)
  % load image
  V = spm_vol(P{pi});
  Y = single(spm_read_vols(V));

  % write unbiased
  Vb = V; [pp,ff,ee] = fileparts(V.fname);
  Vb.fname = fullfile(rdir,sprintf('%s_sd0_amp%0.2f_freq%0.2f%s',ff,0,0,ee)); 
  spm_write_vol(Vb,Y);
  PGT = Vb.fname;

  %% parameter loop - Z
  PTT = cell(1,numel(art.intraSlAmp)*numel(art.interSlFreq)); iPPT = 1;
  for iae = 1:numel(art.intraSlAmp)
    iaa = 1; 
    for ief = 1:numel(art.interSlFreq)
      iaf = 1; %1:numel(art.interSlFreq)

      % add rand distortion
      [Yr,resTr] = cat_vol_resize(Y,'reduceV',1,45,1);
      Yas = 1 - art.intraSlAmp(iaa) + art.intraSlAmp(iaa) * ...
            rand([art.intraSlFreq(iaf) art.intraSlFreq(iaf) size(Y,3)],'single');
      for z=1:size(Yas,3)
         Yas(:,:,z) = Yas(:,:,z) + (1-mean(mean(Yas(:,:,z))));
      end
      resTr.sizeTr(3) = size(Y,3);
      resTr.vx_red(3) = 1;
      Yas = cat_vol_resize(Yas,'dereduceV',resTr);

      % interslicebias
      Za  = 1 - art.interSlAmp(iae) + art.interSlAmp(iae) * ...
            rand([1 1 size(Y,3)],'single');
      Zas = Za - mean(Za(:)); 
      spm_smooth(Zas,Zas,art.interSlFreq(ief)); 
      Zas = Zas .* std(Za(:))./std(Zas(:));
      Zas = Zas + mean(Za(:));
      Zas = repmat(Zas,[size(Y,1),size(Y,2),1]);

      Yb  = Y .* Yas .* Zas;

      % save output
      Vb.fname = fullfile(rdir,sprintf('%s_sdZ_amp%0.2f_freq%0.2f%s',...
                 ff,art.interSlAmp(iae),art.interSlFreq(ief),ee));
      PTT{iPPT} = Vb.fname; iPPT = iPPT + 1;

      spm_write_vol(Vb,Yb);
    end
  end

  %% parameter loop - Y
  PTT = cell(1,numel(art.intraSlAmp)*numel(art.interSlFreq)); iPPT = 1;
  for iae = 1:numel(art.intraSlAmp)
    iaa = 1; 
    for ief = 1:numel(art.interSlFreq)
      iaf = 1; %1:numel(art.interSlFreq)

      % add rand distortion
      [Yr,resTr] = cat_vol_resize(Y,'reduceV',1,45,1);
      Yas = 1 - art.intraSlAmp(iaa) + art.intraSlAmp(iaa) * ...
            rand([art.intraSlFreq(iaf) size(Y,2) art.intraSlFreq(iaf)],'single');
      for y=1:size(Yas,2)
         Yas(:,y,:) = Yas(:,y,:) + (1-mean(mean(Yas(:,y,:))));
      end
      resTr.sizeTr(2) = size(Y,2);
      resTr.vx_red(2) = 1;
      Yas = cat_vol_resize(Yas,'dereduceV',resTr);

      % interslicebias
      Za  = 1 - art.interSlAmp(iae) + art.interSlAmp(iae) * ...
            rand([1 size(Y,2) 1],'single');
      Zas = Za - mean(Za(:)); 
      spm_smooth(Zas,Zas,art.interSlFreq(ief)); 
      Zas = Zas .* std(Za(:))./std(Zas(:));
      Zas = Zas + mean(Za(:));
      Zas = repmat(Zas,[size(Y,1),1,size(Y,3)]);

      Yb  = Y .* Yas .* Zas;

      % save output
      Vb.fname = fullfile(rdir,sprintf('%s_sdY_amp%0.2f_freq%0.2f%s',...
                 ff,art.interSlAmp(iae),art.interSlFreq(ief),ee));
      PTT{iPPT} = Vb.fname; iPPT = iPPT + 1;

      spm_write_vol(Vb,Yb);
    end
  end

  %% parameter loop - X 
  PTT = cell(1,numel(art.intraSlAmp)*numel(art.interSlFreq)); iPPT = 1;
  for iae = 1:numel(art.intraSlAmp)
    iaa = 1; 
    for ief = 1:numel(art.interSlFreq)
      iaf = 1; %1:numel(art.interSlFreq)

      % add rand distortion
      [Yr,resTr] = cat_vol_resize(Y,'reduceV',1,45,1);
      Yas = 1 - art.intraSlAmp(iaa) + art.intraSlAmp(iaa) * ...
            rand([size(Y,1) art.intraSlFreq(iaf) art.intraSlFreq(iaf)],'single');
      for x=1:size(Yas,1)
         Yas(x,:,:) = Yas(x,:,:) + (1-mean(mean(Yas(x,:,:))));
      end
      resTr.sizeTr(1) = size(Y,1);
      resTr.vx_red(1) = 1;
      Yas = cat_vol_resize(Yas,'dereduceV',resTr);

      % interslicebias
      Za  = 1 - art.interSlAmp(iae) + art.interSlAmp(iae) * ...
            rand([size(Y,1) 1 1],'single');
      Zas = Za - mean(Za(:)); 
      spm_smooth(Zas,Zas,art.interSlFreq(ief)); 
      Zas = Zas .* std(Za(:))./std(Zas(:));
      Zas = Zas + mean(Za(:));
      Zas = repmat(Zas,[1,size(Y,2),size(Y,3)]);

      Yb  = Y .* Yas .* Zas;

      % save output
      Vb.fname = fullfile(rdir,sprintf('%s_sdX_amp%0.2f_freq%0.2f%s',...
                 ff,art.interSlAmp(iae),art.interSlFreq(ief),ee));
      PTT{iPPT} = Vb.fname; iPPT = iPPT + 1;

      spm_write_vol(Vb,Yb);
    end
  end
end




%% BWP results
resdir   = '/Users/dahnke/Documents/Dissertation/fig/bwpsliceart';     if ~exist(resdir,'dir'),   mkdir(resdir); end
printdir = '/Users/dahnke/Dropbox/Dissertation/pub/HBM2015/HBM2015_SliceCorrection/bwpsliceart';   if ~exist(printdir,'dir'), mkdir(printdir); end
fig.type = '-depsc'; fig.res  = '-r300';

bwp.dir.Yo  = '/Volumes/vbmDB/MRData/cat12tst/results/deffiles/cat_defaults/BWPC_sliceart_uncorr/';
bwp.dir.Yc  = '/Volumes/vbmDB/MRData/cat12tst/results/deffiles/cat_defaults/BWPC_sliceart_uncorr/';
bwp.dir.Ygt = '/Volumes/MyBook/MRData/cat_tst/+RAW/BWP_Collins/SEG/';
bwp.dir.Yp0o = '/Volumes/vbmDB/MRData/cat12tst/results/deffiles/cat_defaults/BWPC_sliceart_uncorr/';
bwp.dir.Yp0c = '/Volumes/vbmDB/MRData/cat12tst/results/deffiles/cat_defaults/BWPC_sliceart_uncorr/';

Po   = cat_vol_findfiles(bwp.dir.Yo ,'*.nii');
Pc   = cat_vol_findfiles(bwp.dir.Yc ,'scorr*.nii');
Pp0c = cat_vol_findfiles(bwp.dir.Yp0c,'p0scorr*.nii'); 
Pp0o = cat_vol_findfiles(bwp.dir.Yp0o,'p0*.nii'); Pp0o = setdiff(Pp0o,Pp0c);
Pgt  = cat_vol_findfiles(bwp.dir.Ygt,'p0*.nii'); 

% BWP Kappa estimation
[tmp,bwp.K{1}] = cat_tst_calc_kappa(Pp0o,Pgt);
[tmp,bwp.K{2}] = cat_tst_calc_kappa(Pp0c,Pgt);
bwp.Kappa{1}   = mean([bwp.K{1}{2:end-2,3};bwp.K{1}{2:end-2,4}]);
bwp.Kappa{2}   = mean([bwp.K{2}{2:end-2,3};bwp.K{2}{2:end-2,4}]);

% BWP Kappa estimation
bwp.CJV{1} = cat_tst_CJV(Po,Pgt);
bwp.CJV{2} = cat_tst_CJV(Pc,Pgt);

%% analyse
for pi=1:numel(Pp0c)
  [pp,ff,ee] = fileparts(Pp0c{pi});
  bwp.pn(pi)   = str2double(ff(22));
  bwp.rf(pi)   = str2double(ff(27:28));
  bwp.sd(pi)   = ff(48);
  bwp.amp(pi)  = str2double(ff(53:56));
  bwp.freq(pi) = str2double(ff(62:65));
end
%% grouping
if isfield(bwp,'Kgr'), bwp = rmfield(bwp,'Kgr'); bwp = rmfield(bwp,'Knames');end
for gi=1:2
  bwp.Kgr{gi}{1}     = bwp.Kappa{1}(bwp.amp==0.00); bwp.Knames{gi}{1} = 'orig';
  % amp
  bwp.Kgr{gi}{end+1} = bwp.Kappa{gi}(bwp.amp==0.20); bwp.Knames{gi}{end+1} = 'a20';
  bwp.Kgr{gi}{end+1} = bwp.Kappa{gi}(bwp.amp==0.40); bwp.Knames{gi}{end+1} = 'a40';
  bwp.Kgr{gi}{end+1} = bwp.Kappa{gi}(bwp.amp==0.60); bwp.Knames{gi}{end+1} = 'a60';
  % freq
  bwp.Kgr{gi}{end+1} = bwp.Kappa{gi}(bwp.freq==0);   bwp.Knames{gi}{end+1} = 'f0';
  bwp.Kgr{gi}{end+1} = bwp.Kappa{gi}(bwp.freq==2);   bwp.Knames{gi}{end+1} = 'f2';
  bwp.Kgr{gi}{end+1} = bwp.Kappa{gi}(bwp.freq==8);   bwp.Knames{gi}{end+1} = 'f8';
  % dir
  bwp.Kgr{gi}{end+1} = bwp.Kappa{gi}(bwp.sd=='X');   bwp.Knames{gi}{end+1} = 'X';
  bwp.Kgr{gi}{end+1} = bwp.Kappa{gi}(bwp.sd=='Y');   bwp.Knames{gi}{end+1} = 'Y';
  bwp.Kgr{gi}{end+1} = bwp.Kappa{gi}(bwp.sd=='Z');   bwp.Knames{gi}{end+1} = 'Z';
  % qali
  bwp.Kgr{gi}{end+1} = bwp.Kappa{gi}(bwp.pn==3);     bwp.Knames{gi}{end+1} = '3-40';
  bwp.Kgr{gi}{end+1} = bwp.Kappa{gi}(bwp.pn==7);     bwp.Knames{gi}{end+1} = '7-80';
  % corr
  bwp.Kgr{gi}{end+1} = bwp.Kappa{2};                 bwp.Knames{gi}{end+1} = 'corr';
end


%% create groups
if exist('fh1','var') && ishghandle(fh1)
  clf(fh1); figure(fh1)
else
  fh1=figure('Name','figure 1 - bwp main','Position',...
    [0 0 800 400],'color',[1 1 1],'PaperPositionMode','auto');
end

cat_plot_boxplot(bwp.Kgr{1},struct('ylim',[0.3 1],'names',bwp.Knames(1),'title','BWP with slice artifacts',...
  'groupnum',1 ));
%%
print(fh1,fullfile(printdir,sprintf('fig_tst_sliceart_%s',datestr(clock,'yyyymmdd'))),fig.res,fig.type);
%cat_plot_boxplot(bwp.Kgr{2},struct('ylim',[0.85 0.95],'names',bwp.Knames(2),'title','BWP with slice artifacts',...
%  'groupnum',0));

%%
%cat_plot_boxplot(ibsr.CJV,struct('ylim',[0 1])); % not usefull
%cat_plot_boxplot(bwp.Kappa,struct('ylim',[0.5 1]));





%% IBSR results
ibsr.dir.Yo   = '/Volumes/vbmDB/MRData/cat12tst/RAW/IBSRv1';
ibsr.dir.Yc   = '/Volumes/vbmDB/MRData/cat12tst/results/deffiles/cat_defaults/IBSRv1_uncorr/';
ibsr.dir.Ygt  = '/Volumes/MyBook/MRData/cat_tst/+RAW/IBSRv1/SEG/';
ibsr.dir.Yp0o = '/Volumes/vbmDB/MRData/cat12tst/results/deffiles/cat_defaults/IBSRv1_uncorr/';
ibsr.dir.Yp0c = '/Volumes/vbmDB/MRData/cat12tst/results/deffiles/cat_defaults/IBSRv1_uncorr/';

Po   = cat_vol_findfiles(ibsr.dir.Yo ,'*.nii');
Pc   = cat_vol_findfiles(ibsr.dir.Yc ,'scorr*.nii');
Pp0c = cat_vol_findfiles(ibsr.dir.Yp0c,'p0scorr*.nii'); 
Pp0o = cat_vol_findfiles(ibsr.dir.Yp0o,'p0*.nii'); Pp0o = setdiff(Pp0o,Pp0c);
for pi = 1:numel(Pp0o)
  [pp,ff,ee] = fileparts(Pp0o{pi});
  Pgto{pi} = fullfile(ibsr.dir.Ygt,[ff ee]); 
end
for pi = 1:numel(Pp0c)
  [pp,ff,ee] = fileparts(Pp0c{pi});
  Pgtc{pi} = fullfile(ibsr.dir.Ygt,['p0' ff(9:end) ee]); 
end

% BWP Kappa estimation
[tmp,ibsr.K{1}] = cat_tst_calc_kappa(Pp0o,Pgto);
[tmp,ibsr.K{2}] = cat_tst_calc_kappa(Pp0c,Pgtc);
ibsr.Kappa{1}   = mean([ibsr.K{1}{2:end-2,3};ibsr.K{1}{2:end-2,4}]);
ibsr.Kappa{2}   = mean([ibsr.K{2}{2:end-2,3};ibsr.K{2}{2:end-2,4}]);
ibsr.Kappa{2}   = max(mean([ibsr.K{2}{2:end-2,3};ibsr.K{2}{2:end-2,4}]),...
                      mean([ibsr.K{1}{2:end-2,3};ibsr.K{1}{2:end-2,4}])); % use the correction only if necessary

% BWP Kappa estimation
ibsr.CJV{1} = cat_tst_CJV(Po,Pgt);
ibsr.CJV{2} = cat_tst_CJV(Pc,Pgt);

%%
if exist('fh1','var') && ishghandle(fh1)
  clf(fh1); figure(fh1)
else
  fh1=figure('Name','figure 1 - bwp main','Position',...
    [0 0 400 400],'color',[1 1 1],'PaperPositionMode','auto');
end

cat_plot_boxplot(ibsr.Kappa,struct('ylim',[0.7 0.9+eps],'title','BWP with slice artifacts',...
  'groupnum',0,'fontsize',18,'maxwhisker',1));
%%
print(fh1,fullfile(printdir,sprintf('fig_tst_sliceart_ibsr_%s',datestr(clock,'yyyymmdd'))),fig.res,fig.type);
