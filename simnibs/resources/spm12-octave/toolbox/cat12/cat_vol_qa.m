function varargout = cat_vol_qa(action,varargin)
% CAT Preprocessing T1 Quality Assurance
% ______________________________________________________________________
% 
% Estimation of image quality measures like noise, inhomogeneity,
% contrast, resolution, etc. and scaling for school marks. 
%
% [QAS,QAM] = cat_vol_qa(action,varargin)
% 
%
% 1) Use GUI interface to choose segmentation and automatic setting of 
%    original and modified image (if available)
%     [QAS,QAM] = cat_vol_qa()                = cat_vol_qa('p0')
%
%     [QAS,QAM] = cat_vol_qa('p0'[,opt])      - p0 class image
%     [QAS,QAM] = cat_vol_qa('p#'[,opt])      - p1,p2,p3 class images
%     [QAS,QAM] = cat_vol_qa('c#'[,opt])      - c1,c2,c3 class images
%     [QAS,QAM] = cat_vol_qa('*#'[,opt])      - csf,gm,wm class images
%     [QAS,QAM] = cat_vol_qa('p0',Pp0[,opt])           - no GUI call
%     [QAS,QAM] = cat_vol_qa('p#',Pp1,Pp2,Pp3,[,opt])  - no GUI call
%     [QAS,QAM] = cat_vol_qa('c#',Pc1,Pc2,Pc3,[,opt])  - no GUI call
%     [QAS,QAM] = cat_vol_qa('c#',Pcsf,Pgm,Pwm,[,opt]) - no GUI call
%
%
% 2) Use GUI interface to choose all images like for other segmentations
%    and modalities with a similar focus of CSF, GM, and WM tissue 
%    contrast such as PD, T2, or FLASH. 
%     [QAS,QAM] = cat_vol_qa('p0+'[,opt])     - p0 class image  
%     [QAS,QAM] = cat_vol_qa('p#+'[,opt])     - p1,p2,p3 class images  
%     [QAS,QAM] = cat_vol_qa('c#+'[,opt])     - c1,c2,c3 class images 
%     [QAS,QAM] = cat_vol_qa('*#+'[,opt])     - csf,gm,wm class images
%     [QAS,QAM] = cat_vol_qa('p0+',Pp0,Po[,Pm,opt])         - no GUI call
%     [QAS,QAM] = cat_vol_qa('p#+',Pp1,Pp2,Pp3,Po[,Pm,opt]) - no GUI call
%     [QAS,QAM] = cat_vol_qa('c#+',Pc1,Pc2,Pc3,Po[,Pm,opt]) - no GUI call
%
% 
% 3) Use GUI interface to choose all images. I.e. for other segmentations
%    and modalities without focus of GM-WM contrast such as DTI MTI. 
%     [ not implemented yet ]
%
%
% 4) CAT12 internal preprocessing interface 
%    (this is the processing case that is also called in all other cases)
%    [QAS,QAM] = cat_vol_qa('cat12',Yp0,Po,Ym,res[,opt])
%
%
%   Pp0 - segmentation files (p0*.nii)
%   Po  - original files (*.nii)
%   Pm  - modified files (m*.nii)
%   Yp0 - segmentation image matrix
%   Ym  - modified image matrix
%
%   opt            = parameter structure
%   opt.verb       = verbose level  [ 0=nothing | 1=points | 2*=times ]
%   opt.redres     = resolution in mm for intensity scaling [ 4* ];
%   opt.write_csv  = final cms-file
%   opt.write_xml  = images base xml-file
%   opt.sortQATm   = sort QATm output
%     opt.orgval     = original QAM results (no marks)
%     opt.recalc     =
%     opt.avgfactor  = 
%   opt.prefix     = prefix of xml output file (default cat_*.xml) 
%
% ______________________________________________________________________
% Robert Dahnke 
% Structural Brain Mapping Group
% University Jena
%
% $Id: cat_vol_qa.m 1270 2018-02-07 13:32:14Z gaser $
% ______________________________________________________________________

%#ok<*ASGLU>

  % get current release number
  [n, rev_cat] = cat_version;

  % init output
  QAS = struct(); 
  QAR = struct(); 
  cat_qa_warnings = struct('identifier',{},'message',{});
  cat_warnings    = struct('identifier',{},'message',{});
  if nargout>0, varargout = cell(1,nargout); end

  if cat_get_defaults('extopts.subfolders')
    mrifolder = 'mri';
    reportfolder = 'report';
  else
    mrifolder = '';
    reportfolder = '';
  end
   
  % no input and setting of default options
  if nargin==0, action='p0'; end 
  if isstruct(action)
    Pp0 = action.data;
    action = 'p0';
  end
  if nargin>1 && isstruct(varargin{end}) && isstruct(varargin{end})
    opt  = cat_check('checkinopt',varargin{end},defaults);
    nopt = 1; 
  else
    opt  = defaults;
    nopt = 0;
  end

  % check input by action
  switch action
    case {'p0','p0+'}
    % segment image cases
      if nargin<=3
        if (nargin-nopt)<2  
          Pp0 = cellstr(spm_select(inf,'image',...
            'select p0-segment image',{},pwd,'^p0.*')); 
          if isempty(Pp0{1}), return; end
        else
          Pp0 = varargin{1};
        end
        if numel(action)==2
          Po = Pp0; Pm = Pp0;
          for fi=1:numel(Pp0)
            [pp,ff,ee] = spm_fileparts(Pp0{fi});
            [ppa,ppb] = spm_fileparts(pp); 
            if strcmp(ppb,'mri'), ppo = ppa; else ppo = pp; end 

            Po{fi} = fullfile(ppo,[ff(3:end) ee]); 
            Pm{fi} = fullfile(pp,[opt.mprefix  ff(3:end) ee]);
            %Pmv{fi} = fullfile(pp,['m' ff(3:end) ee]); %#ok<AGROW>
            %if ~exist(Pm{fi},'file') && strcmp(opt.mprefix,'nm') && exist(Pmv{fi},'file')
            %  fprintf('Preparing %s.\n',Pmv{fi});
            %  cat_vol_sanlm(Pmv{fi},'n');
            %end

            if ~exist(Po{fi},'file'), Po{fi}=''; end
            if ~exist(Pm{fi},'file'), Pm{fi}=''; end
          end
        else
          Po = cellstr(spm_select(repmat(numel(Pp0),1,2),...
            'image','select original image(s)',{},pwd,'.*')); 
          Pm = cellstr(spm_select(repmat(numel(Pp0),1,2),...
            'image','select modified image(s)',{},pwd,'.*')); 
        end
      elseif nargin<=5
        Pp0 = varargin{1};
        Po  = varargin{2};
        Pm  = varargin{3};
      else
        error('MATLAB:cat_vol_qa:inputerror',...
          'Wrong number/structure of input elements!'); 
      end
    case {'p#','c#','*#','p#+','c#+','*#+'}
    % tissue class image cases
      if nargin-1<=2 % GUI 
        if (nargin-nopt)<2 
          if action(1)=='p' || action(1)=='c'
            % cat/spm case
            Pcsf = cellstr(spm_select(inf,'image',...
              'select p1-segment image',{},pwd,['^' action(1) '1.*'])); 
            if isempty(Pcsf{1}), return; end
            Pgm=Pcsf; Pwm=Pcsf;
            for fi=1:numel(Pcsf)
              [pp,ff,ee] = spm_fileparts(Pcsf{fi});

              Pgm{fi} = fullfile(pp,[action(1) '2' ff(3:end) ee]); 
              Pwm{fi} = fullfile(pp,[action(1) '3' ff(3:end) ee]); 
            end
          else 
            Pcsf = cellstr(spm_select(inf,'image',...
              'select CSF segment image(s)',{},pwd,'.*')); 
            if isempty(Pcsf{1}), return; end
            %Pgm  = cellstr(spm_select(repmat(numel(Pcsf),1,2),...
            %  'image','select GM segment image(s)',{},pwd,'.*')); 
            %Pwm  = cellstr(spm_select(repmat(numel(Pcsf),1,2),...
            %  'image','select WM segment image(s)',{},pwd,'.*')); 
          end 
          if numel(action)==2
            Pp0=Pcsf; Po=Pcsf; Pm=Pcsf;
            for fi=1:numel(Pcsf)
              [pp,ff,ee] = spm_fileparts(Pcsf{fi});
              Po{fi}  = fullfile(pp,[ff(3:end) ee]);
              Pm{fi}  = fullfile(pp,['m'  ff(3:end) ee]);
              Pp0{fi} = fullfile(pp,['p0' ff(3:end) ee]);
            end 
          else
            Po = cellstr(spm_select(repmat(numel(Pcsf),1,2),...
              'image','select original image(s)',{},pwd,'.*')); 
            Pm = cellstr(spm_select(repmat(numel(Pcsf),1,2),...
              'image','select modified image(s)',{},pwd,'.*')); 
            Pp0=Pcsf;
            for fi=1:numel(Pcsf)
              [pp,ff,ee] = spm_fileparts(Pcsf{fi});
              Pp0{fi} = fullfile(pp,['p0' ff(3:end) ee]);
            end 
          end

          % wie komm ich zum p0???
        else
          Pp0 = varargin{1};
        end
      elseif nargin==5 || nargin==6
      else
        error('MATLAB:cat_vol_qa:inputerror',...
          'Wrong number/structure of input elements!'); 
      end

      Yp0 = 1;
    case 'cat12err'
      opt  = cat_check('checkinopt',varargin{end},defaults);
    case 'cat12'
      % CAT12 internal input
      if nargin>3 
        Yp0 = varargin{1};
        Vo  = spm_vol(varargin{2});
        Yo  = single(spm_read_vols(Vo));    
        Ym  = varargin{3}; 
        res = varargin{4};
        V   = res.image;
        cat_warnings = varargin{5};
        species = varargin{6};
        if isfield(varargin{7},'qa')
          if isfield(varargin{7}.qa,'qualitymeasures'), QAR.qualitymeasures = cat_io_updateStruct(QAR,varargin{7}.qa.qualitymeasures); end
          if isfield(varargin{7}.qa,'subjectmeasures'), QAS.subjectmeasures = cat_io_updateStruct(QAS,varargin{7}.qa.subjectmeasures); end
        end
        % opt = varargin{end} in line 96)
        opt.verb = 0;
        
        % reduce to original native space if it was interpolated
        if any(size(Yp0)~=Vo.dim)
          if isfield(Vo,'private'), Vo = rmfield(Vo,'private'); end
          if isfield(Vo,'mat0'),    Vo = rmfield(Vo,'mat0');    end
          Vo.dat = zeros(Vo.dim,'single'); Vo.dt(1) = 16; Vo.pinfo = [1;0;0];
          
          Vp0t          = res.image;
          if isfield(Vp0t,'private'), Vp0t = rmfield(Vp0t,'private'); end
          if isfield(Vp0t,'mat0'),    Vp0t = rmfield(Vp0t,'mat0'); end
          Vp0t.dt(1)    = 16;
          Vp0t.pinfo    = [1;0;0];
          Vp0t.dat      = Yp0;

          % resampling and corrections of the Yp0
         % Vp0t       = spm_write_vol(Vp0t,double(Yp0));
          [Vtpm,Yp0] = cat_vol_imcalc(Vp0t,Vo,'i1',struct('interp',2,'verb',0));
          rf         = 50;
          Yp0        = single(Yp0);
          Yp0r       = round(Yp0*rf)/rf;
          YMR        = false(size(Yp0));
          for i=1:4, YMR = YMR | (Yp0>(i-1/rf) & Yp0<(i+1/rf)); end
          Yp0(YMR)   = Yp0r(YMR); clear YMR Ynr;
          
          % resampling of the corrected image
          Vp0t.dat   = Ym;
          [Vtpm,Ym]  = cat_vol_imcalc(Vp0t,Vo,'i1',struct('interp',6,'verb',0)); 
          Ym         = single(Ym);
        end
        
      else
        error('MATLAB:cat_vol_qa:inputerror',...
          'Wrong number/structure of input elements!'); 
      end
    otherwise
      error('MATLAB:cat_vol_qa:inputerror',...
        'Wrong number/structure of input elements!'); 
  end
  if ~exist('species','var'), species='human'; end
    
  
  %
  % --------------------------------------------------------------------
  [QA,QMAfn]  = cat_stat_marks('init'); 
  stime  = clock;
  
  
  
  % Print options
  % --------------------------------------------------------------------
  opt.snspace = [70,7,3];
  Cheader = {'scan'};
  Theader = sprintf(sprintf('%%%ds:',opt.snspace(1)-1),'scan');
  Tline   = sprintf('%%5d) %%%ds:',opt.snspace(1)-8);
  Tline2  = sprintf('%%5d) %%6s%%%ds:',opt.snspace(1)-14); 
  Tavg    = sprintf('%%%ds:',opt.snspace(1)-1);
  TlineE  = sprintf('%%5d) %%%ds: %%s',opt.snspace(1)-7);
  for fi=1:numel(QMAfn)
    Cheader = [Cheader QMAfn{fi}]; %#ok<AGROW>
    Theader = sprintf(sprintf('%%s%%%ds',opt.snspace(2)),Theader,...
                QMAfn{fi}(1:min(opt.snspace(2)-1,numel(QMAfn{fi}))));
    Tline   = sprintf('%s%%%d.%df',Tline,opt.snspace(2),opt.snspace(3));
    Tline2  = sprintf('%s%%%d.%df',Tline2,opt.snspace(2),opt.snspace(3));
    Tavg    = sprintf('%s%%%d.%df',Tavg,opt.snspace(2),opt.snspace(3));
  end
  Cheader = [Cheader 'IQM'];
  Theader = sprintf(sprintf('%%s%%%ds',opt.snspace(2)),Theader,'IQM');
  Tline   = sprintf('%s%%%d.%df\n',Tline,opt.snspace(2),opt.snspace(3));
  Tline2  = sprintf('%s%%%d.%df\n',Tline2,opt.snspace(2),opt.snspace(3));
  Tavg    = sprintf('%s%%%d.%df\n',Tavg,opt.snspace(2),opt.snspace(3));
  
  
  

  
  
  % estimation part    
  switch action
    case {'p0','p#','c#','*#','p0+','p#+','c#+','*#+'}    
    % loop for multiple files 
      % return for empty input
      if isempty(Pp0) || (isempty(Pp0{1}) && numel(Pp0)<=1) 
        cat_io_cprintf('com','No images for QA!\n'); 
        return
      end
      
      if opt.verb>1
        fprintf('\n%s\n\n%s\n%s\n', ...
          sprintf('CAT Preprocessing T1 Quality Assurance (%s):',...
          sprintf('Rev: %s',rev_cat)), Theader,repmat('-',size(Theader)));  
      end

      qamat   = nan(numel(Po),numel(QMAfn));
      qamatm  = nan(numel(Po),numel(QMAfn));
      mqamatm = 10.5*ones(numel(Po),1);
    
      
      QAS = struct(); QAR = struct(); 
      QAR.mark2rps = @(mark) min(100,max(0,105 - mark*10));
      
      for fi=1:numel(Pp0)
        try
          if exist(Po{fi},'file')
            Vo  = spm_vol(Po{fi});
          else
            error('cat_vol_qa:noYo','No original image.');
          end
  
          Yp0 = single(spm_read_vols(spm_vol(Pp0{fi})));
          if ~isempty(Pm{fi}) && exist(Pm{fi},'file')
            Ym  = single(spm_read_vols(spm_vol(Pm{fi})));
          else
            error('cat_vol_qa:noYm','No corrected image.');
          end
  
          res.image = spm_vol(Pp0{fi}); 
          [QASfi,QAMfi,cat_qa_warnings{fi}] = cat_vol_qa('cat12',Yp0,Vo,Ym,res,cat_warnings,species,opt);

     
          QAS = cat_io_updateStruct(QAS,QASfi,0,fi);
          QAR = cat_io_updateStruct(QAR,QAMfi,0,fi);
        
          
          % color for the differen mark cases (opt.process)
          for fni=1:numel(QMAfn)
            qamat(fi,fni)  = QAS(fi).qualitymeasures.(QMAfn{fni});
            qamatm(fi,fni) = QAR(fi).qualityratings.(QMAfn{fni});
          end
          mqamatm(fi) = QAR(fi).qualityratings.IQR;
          mqamatm(fi) = max(0,min(10.5, mqamatm(fi)));
          
          
          % print the results for each scan 
          if opt.verb>1 
            if opt.orgval 
              cat_io_cprintf(opt.MarkColor(max(1,round( mqamatm(fi,:)/9.5 * ...
                size(opt.MarkColor,1))),:),sprintf(Tline,fi,...
                QAS(fi).filedata.fnames, ... spm_str_manip(QAS(fi).filedata.file,['f' num2str(opt.snspace(1) - 14)]),...
                qamat(fi,:),max(1,min(6,mqamatm(fi)))));
            else
              cat_io_cprintf(opt.MarkColor(max(1,round( mqamatm(fi,:)/9.5 * ...
                size(opt.MarkColor,1))),:),sprintf(Tline,fi,...
                QAS(fi).filedata.fnames, ... spm_str_manip(QAS(fi).filedata.file,['f' num2str(opt.snspace(1) - 14)]),...
                qamatm(fi,:),max(1,min(6,mqamatm(fi)))));
            end
          end
        catch  %#ok<CTCH> ... normal "catch err" does not work for MATLAB 2007a
          try %#ok<TRYNC>
            e = lasterror; %#ok<LERR> ... normal "catch err" does not work for MATLAB 2007a
         
            switch e.identifier
              case {'cat_vol_qa:noYo','cat_vol_qa:noYm','cat_vol_qa:badSegmentation'}
                em = e.identifier;
              otherwise
                em = ['ERROR:\n' repmat(' ',1,10) e.message '\n'];
                for ei=1:numel(e.stack)
                  em = sprintf('%s%s%5d: %s\n',em,repmat(' ',1,10),...
                    e.stack(ei).line(end),e.stack(ei).name);
                end  
            end

            [pp,ff] = spm_fileparts(Po{fi});
            QAS(fi).filedata.fnames = [spm_str_manip(pp,sprintf('k%d',floor( (opt.snspace(1)-19) /3) - 1)),'/',...
                                 spm_str_manip(ff,sprintf('k%d',(opt.snspace(1)-19) - floor((opt.snspace(1)-14)/3)))];
            cat_io_cprintf(opt.MarkColor(end,:),sprintf(TlineE,fi,...
               QAS(fi).filedata.fnames,[em '\n']));
%            spm_str_manip(Po{fi},['f' num2str(opt.snspace(1) - 14)]),em));
          end
        end
      end      
      
      
     
      % sort by mean mark
      % ----------------------------------------------------------------
      if opt.sortQATm && numel(Po)>1
        % sort matrix
        [smqamatm,smqamatmi] = sort(mqamatm,'ascend');
        sqamatm  = qamatm(smqamatmi,:);
        sqamat   = qamat(smqamatmi,:); 

        % print matrix
        if opt.verb>0
          fprintf('%s\n',repmat('-',size(Theader))); 
          for fi=1:numel(QAS)
            if opt.orgval 
              cat_io_cprintf(opt.MarkColor(min(size(opt.MarkColor,1),...
                round( smqamatm(fi,:)/9.5 * ...
                size(opt.MarkColor,1))),:),sprintf(...
                Tline2,fi,sprintf('(%d)',smqamatmi(fi)),...
                QAS(smqamatmi(fi)).filedata.fnames, ...
                ...spm_str_manip(QAS(smqamatmi(fi)).filedata.file,['f' num2str(opt.snspace(1) - 14)]),...
                sqamat(fi,:),max(1,min(6,smqamatm(fi)))));
            else
              cat_io_cprintf(opt.MarkColor(min(size(opt.MarkColor,1),...
                round( smqamatm(fi,:)/9.5 * ...
                size(opt.MarkColor,1))),:),sprintf(...
                Tline2,fi,sprintf('(%d)',smqamatmi(fi)),...
                QAS(smqamatmi(fi)).filedata.fnames, ...
                ...spm_str_manip(QAS(smqamatmi(fi)).filedata.file,['f' num2str(opt.snspace(1) - 14)]),...
                sqamatm(fi,:),smqamatm(fi)));
            end
          end
        end
      else
        %[smqamatm,smqamatmi] = sort(mqamatm,'ascend');
        %sqamatm  = qamatm(smqamatmi,:);
      end
      % print the results for each scan 
      if opt.verb>1 && numel(Pp0)>1
        fprintf('%s\n',repmat('-',size(Theader)));  
        if opt.orgval 
          fprintf(Tavg,'mean',cat_stat_nanmean(qamat,1),cat_stat_nanmean(mqamatm,1));    %#ok<CTPCT>
          fprintf(Tavg,'std' , cat_stat_nanstd(qamat,1), cat_stat_nanstd(mqamatm,1));    %#ok<CTPCT>  
        else
          fprintf(Tavg,'mean',cat_stat_nanmean(qamatm,1),cat_stat_nanmean(mqamatm,1));    %#ok<CTPCT>
          fprintf(Tavg,'std' , cat_stat_nanstd(qamatm,1), cat_stat_nanstd(mqamatm,1));    %#ok<CTPCT>  
        end 
        %fprintf('%s\n',repmat('-',size(Theader)));  
        %fprintf(Tavg,'mean',mean(qamat,1));  
        %fprintf(Tavg,'std', std(qamat,1));    
      end
      if opt.verb>0, fprintf('\n'); end


      
      % result tables (cell structures)
      % ----------------------------------------------------------------
      if nargout>2 && opt.write_csv
        QAT   = [Cheader(1:end-1); ... there is no mean for the original measures
                 Po               , num2cell(qamat); ...
                 'mean'           , num2cell(cat_stat_nanmean(qamat,1)); ...
                 'std'            , num2cell( cat_stat_nanstd(qamat,1,1))];
        QATm  = [Cheader; ...
                 Po               , num2cell(qamatm)          , ...
                                    num2cell(cat_stat_nanmean(qamatm,2)); ...
                 'mean'           , num2cell(cat_stat_nanmean(qamatm,1))  , ...
                                    num2cell(cat_stat_nanmean(mqamatm,1)); ...
                 'std'            , num2cell( cat_stat_nanstd(qamatm,1,1)), ...
                                    num2cell( cat_stat_nanstd(mqamatm,1))];


        % write csv results
        % --------------------------------------------------------------
        if opt.write_csv
          pp = spm_fileparts(Pp0{1});
          cat_io_csv(fullfile(pp,reportfolder,[opt.prefix num2str(numel(Vo),'%04d') ...
            'cat_vol_qa_values.csv']),QAT);
          cat_io_csv(fullfile(pp,reportfolder,[opt.prefix num2str(numel(Vo),'%04d') ...
            'cat_vol_qa_marks.csv']),QATm);
        end
      end 
      
      if opt.verb>0
        fprintf('Quality Control for %d subject was done in %0.0fs\n', ...
          numel(Pp0),etime(clock,stime)); fprintf('\n');
      end
  
      
    case 'cat12err'
      
      % file information
      % ----------------------------------------------------------------
      [pp,ff,ee] = spm_fileparts(opt.job.channel.vols{opt.subj});
      [QAS.filedata.path,QAS.filedata.file] = spm_fileparts(opt.job.channel.vols{opt.subj});
      QAS.filedata.fname  = opt.job.data{opt.subj};
      QAS.filedata.F      = opt.job.data{opt.subj}; 
      QAS.filedata.Fm     = fullfile(pp,mrifolder,['m'  ff ee]);
      QAS.filedata.Fp0    = fullfile(pp,mrifolder,['p0' ff ee]);
      QAS.filedata.fnames = [spm_str_manip(pp,sprintf('k%d',...
                         floor( max(opt.snspace(1)-19-ff,opt.snspace(1)-19)/3) - 1)),'/',...
                       spm_str_manip(ff,sprintf('k%d',...
                         (opt.snspace(1)-19) - floor((opt.snspace(1)-14)/3)))];
    

      % software, parameter and job information
      % ----------------------------------------------------------------
      [nam,rev_spm] = spm('Ver');
      QAS.software.version_spm = rev_spm;
      A = ver;
      for i=1:length(A)
        if strcmp(A(i).Name,'MATLAB'),
          QAS.software.version_matlab = A(i).Version; 
        end
      end
      clear A
      QAS.software.version_cat  = rev_cat;
      QAS.software.function     = which('cat_vol_qa');
      QAS.software.markdefs     = which('cat_stat_marks');
      QAS.software.qamethod     = action; 
      QAS.software.date         = datestr(clock,'yyyymmdd-HHMMSS');
      warning off
      QAS.software.opengl       = opengl('INFO');
      QAS.software.opengldata   = opengl('DATA');
      warning on
      QAS.hardware.computer     = mexext; 
      try
        QAS.hardware.numcores = max(feature('numcores'),1);
      catch
        QAS.hardware.numcores = 1;
      end
      
      
      %opt.job  = rmfield(opt.job,{'data','channel','output'}); 
      %QAS.parameter             = opt.job; 
      QAS.parameter.opts        = opt.job.opts;
      QAS.parameter.extopts     = opt.job.extopts;
      %QAS.parameter.output      = opt.job.output;
      QAS.parameter.caterr      = opt.caterr; 
      QAS.error                 = opt.caterrtxt; 
      
      % export 
      if opt.write_xml
        cat_io_xml(fullfile(pp,reportfolder,[opt.prefix ff '.xml']),QAS,'write');
      end
      
    case 'cat12'
    % estimation of the measures for the single case    
    
    
      % file information
      % ----------------------------------------------------------------
      [pp,ff,ee] = spm_fileparts(Vo.fname);
      [QAS.filedata.path,QAS.filedata.file] = spm_fileparts(Vo.fname);
      QAS.filedata.fname  = Vo.fname;
      QAS.filedata.F      = Vo.fname; 
      QAS.filedata.Fm     = fullfile(pp,mrifolder,['m'  ff ee]);
      QAS.filedata.Fp0    = fullfile(pp,mrifolder,['p0' ff ee]);
      QAS.filedata.fnames = [spm_str_manip(pp,sprintf('k%d',...
                         floor( max(opt.snspace(1)-19-ff,opt.snspace(1)-19)/3) - 1)),'/',...
                       spm_str_manip(ff,sprintf('k%d',...
                         (opt.snspace(1)-19) - floor((opt.snspace(1)-14)/3)))];
    

      % software, parameter and job information
      % ----------------------------------------------------------------
      [nam,rev_spm] = spm('Ver');
      QAS.software.version_spm = rev_spm;
      A = ver;
      for i=1:length(A)
        if strcmp(A(i).Name,'MATLAB'),
          QAS.software.version_matlab = A(i).Version; 
        end
      end
      clear A
      QAS.software.version_cat  = rev_cat;
      QAS.software.function     = which('cat_vol_qa');
      QAS.software.markdefs     = which('cat_stat_marks');
      QAS.software.qamethod     = action; 
      QAS.software.date         = datestr(clock,'yyyymmdd-HHMMSS');
      warning off
      QAS.software.opengl       = opengl('INFO');
      QAS.software.opengldata   = opengl('DATA');
      warning on
      QAS.software.cat_warnings = cat_warnings;
 
      %QAS.parameter             = opt.job; 
      if isfield(opt,'job')
        QAS.parameter.opts        = opt.job.opts;
        QAS.parameter.extopts     = opt.job.extopts;
        %QAS.parameter.output      = opt.job.output;
        if exist('res','var');
          rf = {'Affine','lkp','mn','vr'}; % important SPM preprocessing variables
          for rfi=1:numel(rf)
            if isfield(res,rf{rfi}), QAS.parameter.spm.(rf{rfi}) = res.(rf{rfi}); end
          end
        end
      end
     
      %% resolution, boundary box
      %  ---------------------------------------------------------------
      QAS.software.cat_qa_warnings = struct('identifier',{},'message',{});
      vx_vol  = sqrt(sum(Vo.mat(1:3,1:3).^2));
      vx_voli = sqrt(sum(V.mat(1:3,1:3).^2));
      Yp0toC  = @(Yp0,c) 1-min(1,abs(Yp0-c));
      
      %  resolution 
      QAS.qualitymeasures.res_vx_vol    = vx_vol;
      if 1 % CAT internal resolution
        QAS.qualitymeasures.res_vx_voli = vx_voli;
      end
      QAS.qualitymeasures.res_RMS       = cat_stat_nanmean(vx_vol.^2).^0.5;
      % further unused measure (just for test/comparison)
      %QAS.qualitymeasures.res_isotropy  = max(vx_vol)./min(vx_vol);
      %QAS.qualitymeasures.res_vol       = prod(abs(vx_vol));
      %QAS.qualitymeasures.res_MVR       = mean(vx_vol);
      
      % boundary box - brain tissue next to image boundary
      bbth = round(2/cat_stat_nanmean(vx_vol)); M = true(size(Yp0));
      M(bbth:end-bbth,bbth:end-bbth,bbth:end-bbth) = 0;
      QAS.qualitymeasures.res_BB = sum(Yp0(:)>1.25 & M(:))*prod(abs(vx_vol)); 

      % check segmentation
      spec = species; for ai=num2str(0:9); spec = strrep(spec,ai,''); end; 
      bvol = species; for ai=char(65:122); bvol = strrep(bvol,ai,''); end; bvol = str2double(bvol);
      
      subvol = [sum(Yp0(:)>2.5 & Yp0(:)<3.1)*prod(vx_vol)/1000,... 
                sum(Yp0(:)>1.5 & Yp0(:)<2.5)*prod(vx_vol)/1000,...
                sum(Yp0(:)>0.5 & Yp0(:)<1.5)*prod(vx_vol)/1000]; 
      
      if isempty(bvol) 
        switch spec
          case 'human'
            bvol = 1400; 
          otherwise
            warning('cat_vol_qa:species',...
              sprintf('Unknown species %s (C=%0.0f,G=%0.0f,W=%0.0f).',species,subvol)); %#ok<SPWRN>
        end
      end
      if  sum(subvol)<bvol/3 || sum(subvol)>bvol*3
        warning('cat_vol_qa:badSegmentation',...
          sprintf('Bad %s segmentation (C=%0.0f,G=%0.0f,W=%0.0f).',species,subvol)) %#ok<SPWRN>
      end
      

      %%  estimate QA
      %  ---------------------------------------------------------------
      % remove space arount the brain for speed-up
      [Yo,Ym,Yp0]   = cat_vol_resize({Yo,Ym,Yp0},'reduceBrain',vx_vol,4,Yp0>1.5);
      
      % rought contast and noise estimation to get a stable T1 map for threshold estimation
      T1th = [cat_stat_nanmedian(Ym(Yp0toC(Yp0(:),1)>0.9)) ...
              cat_stat_nanmedian(Ym(Yp0toC(Yp0(:),2)>0.9)) ...
              cat_stat_nanmedian(Ym(Yp0toC(Yp0(:),3)>0.9))];
      noise = max(0,min(1,cat_stat_nanstd(Ym(Yp0(:)>2.9)) / min(abs(diff(T1th)))));
      Yms = Ym+0; spm_smooth(Yms,Yms,repmat(double(noise)*4,1,3));      % smoothing to reduce high frequency noise
            
      % basic tissue classes - erosion to avoid PVE, std to avoid other tissues (like WMHs)
      voli = @(v) (v ./ (pi * 4./3)).^(1/3); 
      rad  = voli( QAS.subjectmeasures.vol_TIV) ./ cat_stat_nanmean(vx_vol);
      Ysc  = 1-cat_vol_smooth3X(Yp0<1 | Yo==0,min(24,max(16,rad*2)));   % fast 'distance' map
      Ycm  = cat_vol_morph(Yp0>0.5 & Yp0<1.5 & Yms<cat_stat_nanmean(T1th(1:2)),'e') & ...
              Ysc>0.75 & Yp0<1.25;% avoid PVE & ventricle focus
      if sum(Ycm(:)>0)<10; Ycm=cat_vol_morph(Yp0>0.5 & Yp0<1.5 & Yms<cat_stat_nanmean(T1th(1:2)),'e') & Yp0<1.25; end
      if sum(Ycm(:)>0)<10; Ycm=Yp0>0.5 & Yms<cat_stat_nanmean(T1th(1:2)) & Yp0<1.25; end
      %Ycm  = Ycm | (Yp0==1 & Ysc>0.7 & Yms<cat_stat_nanmean(T1th(2:3))); % HEBEL      
      Ygm1 = round(Yp0*10)/10==2;                                       % avoid PVE 1
      Ygm2 = cat_vol_morph(Yp0>1.1,'e') & cat_vol_morph(Yp0<2.9,'e');   % avoid PVE 2
      Ygm  = (Ygm1 | Ygm2) & Ysc<0.9;                                   % avoid PVE & no subcortex
      Ywm  = cat_vol_morph(Yp0>2.1,'e') & Yp0>2.9 & ...                 % avoid PVE & subcortex
        Yms>min(cat_stat_nanmean(T1th(2:3)),(T1th(2) + 2*noise*diff(T1th(2:3))));   % avoid WMHs2
      clear Ygm1 Ygm2; % Ysc; 
      
      %% further refinements of the tissue maps
      T2th = [median(Yms(Ycm)) median(Yms(Ygm)) median(Yms(Ywm))];
      Ycm  = Ycm & Yms>(T2th(1)-16*noise*diff(T2th(1:2))) & Ysc &...
             Yms<(T2th(1)+0.1*noise*diff(T2th(1:2)));
      if sum(Ycm(:)>0)<10; Ycm=cat_vol_morph(Yp0>0.5 & Yp0<1.5 & Yms<cat_stat_nanmean(T1th(1:2)),'e') & Yp0<1.25; end
      if sum(Ycm(:)>0)<10; Ycm=Yp0>0.5 & Yms<cat_stat_nanmean(T1th(1:2)) & Yp0<1.25; end     
      Ygm  = Ygm & Yms>(T2th(2)-2*noise*diff(T1th(2:3))) & Yms<(T2th(2)+2*noise*diff(T1th(2:3)));
      Ygm(smooth3(Ygm)<0.2) = 0;
      Ycm  = cat_vol_morph(Ycm,'lc'); % to avoid holes
      Ywm  = cat_vol_morph(Ywm,'lc'); % to avoid holes
      Ywe  = cat_vol_morph(Ywm,'e');  
      
      
      %% low resolution tissue intensity maps (smoothing)
      % High frequency noise is mostly uncritical as far as simple smoothing can reduce it. 
      % Although the very low frequency interferences (inhomogeneity) is unproblematic in most cases,  
      % but will influence the noise pattern. 
      % But most important is the noise with the medium high frequencies, that we try do detect by 
      % reducing the very high and low noise pattern by filtering and pixel smoothing by reduction.
      res = 2; vx_volx = 1; 
      Yos = cat_vol_localstat(Yo,Ywm,1,1); Yo(Yos>0)=Yos(Yos>0);        % reduce high frequency noise in WM 
      Yos = cat_vol_localstat(Yo,Ycm,1,1); Yo(Yos>0)=Yos(Yos>0);        % reduce high frequency noise in CSF

      Yc  = cat_vol_resize(Yo .* Ycm,'reduceV',vx_volx,res,32,'min');   % CSF thr. (minimum to avoid PVE)
      Yg  = cat_vol_resize(Yo .* Ygm,'reduceV',vx_volx,res,32,'meanm'); % GM thr.
      Yw  = cat_vol_resize(Yo .* Ywe,'reduceV',vx_volx,res,32,'meanm'); % WM thr. and bias correction (Ywme)
      Ywc = cat_vol_resize(Ym .* Ywe,'reduceV',vx_volx,res,32,'meanm'); % for bias correction
      Ywb = cat_vol_resize(Yo .* Ywm,'reduceV',vx_volx,res,32,'max');   % for WM inhomogeneity estimation (avoid PVE)
      Ywn = cat_vol_resize(Yo .* Ywm,'reduceV',vx_volx,res,32,'meanm'); % for WM noise
      Ycn = cat_vol_resize(Yo .* Ycm,'reduceV',vx_volx,res,32,'meanm'); % for CSF noise
      Ycm = cat_vol_resize(Ycm      ,'reduceV',vx_volx,res,32,'meanm'); % CSF thr. (minimum to avoid PVE)
      Ygm = cat_vol_resize(Ygm      ,'reduceV',vx_volx,res,32,'meanm'); % GM thr.
      Ywm = cat_vol_resize(Ywm      ,'reduceV',vx_volx,res,32,'meanm'); % WM thr. and bias correction (Ywme)
      Ywe = cat_vol_resize(Ywe      ,'reduceV',vx_volx,res,32,'meanm'); % WM thr. and bias correction (Ywme)
      
      % only voxel that were the product of 
      Yc  = Yc  .* (Ycm>=0.5); Yg  = Yg  .* (Ygm>=0.5);  Yw  = Yw  .* (Ywe>=0.5); 
      Ywc = Ywc .* (Ywe>=0.5); Ywb = Ywb .* (Ywm>=0.5);  Ywn = Ywn .* (Ywm>=0.5);
      Ycn = Ycn .* (Ycm>=0.5);
      
      clear Ycm Ygm Ywm Ywme;
      [Yo,Ym,Yp0,resr] = cat_vol_resize({Yo,Ym,Yp0},'reduceV',vx_volx,res,32,'meanm'); 
      resr.vx_volo = vx_vol; vx_vol=resr.vx_red .* resr.vx_volo;
      
      % intensity scaling for normalized Ym maps like in CAT12
      Ywc = Ywc .* (cat_stat_nanmean(Yo(Yp0(:)>2))/cat_stat_nanmean(Ym(Yp0(:)>2)));
      
      %% bias correction for original map, based on the 
      WI  = Yw./max(eps,Ywc); WI(isnan(WI) | isinf(WI)) = 0; 
      WI  = cat_vol_approx(WI,2);
      WI  = cat_vol_smooth3X(WI,1);
      Ywn = Ywn./WI; Ywn = round(Ywn*1000)/1000;
      Ymi = Yo ./WI; Ymi = round(Ymi*1000)/1000;
      Yc  = Yc ./WI; Yc  = round(Yc *1000)/1000;
      Yg  = Yg ./WI; Yg  = round(Yg *1000)/1000;
      Yw  = Yw ./WI; Yw  = round(Yw *1000)/1000;
      clear WIs ;
      
      Ywb = Ywb ./ cat_stat_nanmean(Ywb(Yp0(:)>2));
     
      % tissue segments for contrast estimation etc. 
      CSFth = cat_stat_nanmean(Yc(~isnan(Yc(:)) & Yc(:)~=0)); 
      GMth  = cat_stat_nanmean(Yg(~isnan(Yg(:)) & Yg(:)~=0));
      WMth  = cat_stat_nanmean(Yw(~isnan(Yw(:)) & Yw(:)~=0)); 
      T3th  = [CSFth GMth WMth];
      
      % estimate background
      [Ymir,resYbg] = cat_vol_resize(Ymi,'reduceV',1,6,32,'meanm'); 
      try
        warning 'off' 'MATLAB:cat_vol_morph:NoObject'
        BGCth = min(T3th)/2; 
        Ybgr = cat_vol_morph(cat_vol_morph(Ymir<BGCth,'lc',1),'e',2/cat_stat_nanmean(resYbg.vx_volr)) & ~isnan(Ymir);
        Ybg  = cat_vol_resize(Ybgr,'dereduceV',resYbg)>0.5; clear Yosr Ybgr;
        if sum(Ybg(:))<32, Ybg = cat_vol_morph(Yo<BGCth,'lc',1) & ~isnan(Yo); end
        warning 'on'  'MATLAB:cat_vol_morph:NoObject'
        BGth = cat_stat_nanmedian(Ymi(Ybg(:)));   
      catch
        warning 'on'  'MATLAB:cat_vol_morph:NoObject'
        try
          % non-zero background
          Ygr  = cat_vol_grad(Ymir); 
          warning 'off' 'MATLAB:cat_vol_morph:NoObject'
          Ybgr = cat_vol_morph(cat_vol_morph(Ygr<0.3 & Yp0<0,'lc',1),'e',2/cat_stat_nanmean(resYbg.vx_volr)) & ~isnan(Ymir);
          Ybg  = cat_vol_resize(Ybgr,'dereduceV',resYbg)>0.5; clear Yosr Ybgr;
          if sum(Ybg(:))<32, Ybg = cat_vol_morph(Yo<BGCth,'lc',1) & ~isnan(Yo); end
          warning 'on'  'MATLAB:cat_vol_morph:NoObject'
          BGth = cat_stat_nanmedian(Ymi(Ybg(:)));   
        catch
          warning 'on'  'MATLAB:cat_vol_morph:NoObject'
          BGth = nan; 
        end
      end
          
      % (relative) average tissue intensity of each class
      QAS.qualitymeasures.tissue_mn  = ([BGth CSFth GMth WMth]);
      QAS.qualitymeasures.tissue_mnr = QAS.qualitymeasures.tissue_mn ./ (max([WMth,GMth])); 
      if WMth>GMth
        QAS.qualitymeasures.tissue_weighting = 'T1';
      elseif WMth<GMth && GMth<CSFth
        QAS.qualitymeasures.tissue_weighting = 'inverse';
      end
      
      % (relative) standard deviation of each class
      QAS.qualitymeasures.tissue_std(1) = cat_stat_nanstd( Ymi(Ybg(:)) );
      for ci=2:4
        QAS.qualitymeasures.tissue_std(ci) = cat_stat_nanstd(Ymi(Yp0toC(Yp0(:),ci-1)>0.5 & ~isinf(Yp0(:))));
      end
      QAS.qualitymeasures.tissue_stdr = QAS.qualitymeasures.tissue_std ./ (WMth-BGth);
       
      % (relative) (mininum) tissue contrast ( CSF-GM-WM ) 
      % - the CSF threshold varies strongly due to bad segmentations,
      %   and anatomica variance, so its better to use GM-WM contrast 
      %   and take care of overoptimisation with values strongly >1/3
      %   of the relative contrast
      contrast = min(abs(diff(QAS.qualitymeasures.tissue_mn(3:4)))) ./ (max([WMth,GMth])); % default contrast
      contrast = contrast + min(0,13/36 - contrast)*1.2;                      % avoid overoptimsization
      QAS.qualitymeasures.contrast  = contrast * (max([WMth,GMth])); 
      QAS.qualitymeasures.contrastr = contrast;

      
      
      %% noise estimation (original (bias corrected) image)
      % WM variance only in one direction to avoid WMHs!
      rms=1; nb=1;
      NCww = sum(Ywn(:)>0) * prod(vx_vol);
      NCwc = sum(Ycn(:)>0) * prod(vx_vol);
      [Yos2,YM2] = cat_vol_resize({Ywn,Ywn>0},'reduceV',vx_vol,3,16,'meanm');
      NCRw = estimateNoiseLevel(Yos2,YM2>0.5,nb,rms) / max(GMth,WMth) / contrast  ; 
      if BGth<-0.1 && WMth<3, NCRw=NCRw/3; end% MT weighting
      clear Yos0 Yos1 Yos2 YM0 YM1 YM2;
        
      % CSF variance of large ventricle
      % for typical T2 images we have to much signal in the CSF and can't use it for noise estimation!
      wcth = 200; 
      if CSFth<GMth && NCwc>wcth
        [Yos2,YM2] = cat_vol_resize({Ycn,Ycn>0},'reduceV',vx_vol,3,16,'meanm');
        NCRc = estimateNoiseLevel(Yos2,YM2>0.5,nb,rms) / max(GMth,WMth)  / contrast ; 
        clear Yos0 Yos1 Yos2 YM0 YM1 YM2;
      else
        NCRc = 0;
        NCwc = 0;
      end
      % 1/sqrt(volume) to compensate for noise differency due to different volumen size. 
      % Overall there are better chances to correct high resolution noise.
      % Nitz W R. Praxiskurs MRT. Page 28.
      NCwc = min(wcth,max(0,NCwc-wcth)); NCww = min(wcth,NCww) - NCwc; % use CSF if possible
      if NCwc<3*wcth && NCww<10*wcth, NCRc = min(NCRc,NCRw); end
      QAS.qualitymeasures.NCR = (NCRw*NCww + NCRc*NCwc)/(NCww+NCwc);
      QAS.qualitymeasures.NCR = QAS.qualitymeasures.NCR * (prod(resr.vx_volo*res))^0.4 * 5/4; %* 7.5; %15;
      %QAS.qualitymeasures.CNR = 1 / QAS.qualitymeasures.NCR;  
%fprintf('NCRw: %8.3f, NCRc: %8.3f, NCRf: %8.3f\n',NCRw,NCRc,(NCRw*NCww + NCRc*NCwc)/(NCww+NCwc));

      
      %% Bias/Inhomogeneity (original image with smoothed WM segment)
      Yosm = cat_vol_resize(Ywb,'reduceV',vx_vol,3,32,'meanm');      % resolution and noise reduction
      for si=1:max(1,min(3,round(QAS.qualitymeasures.NCR*4))), Yosm = cat_vol_localstat(Yosm,Yosm>0,1,1); end 
      QAS.qualitymeasures.ICR  = cat_stat_nanstd(Yosm(Yosm(:)>0)) / contrast;
      %QAS.qualitymeasures.CIR  = 1 / QAS.qualitymeasures.ICR;

      
  
      %% marks
      QAR = cat_stat_marks('eval',1,QAS);

      % export 
      if opt.write_xml
        QAS.qualityratings = QAR.qualityratings;
        QAS.subjectratings = QAR.subjectratings;
        
        cat_io_xml(fullfile(pp,reportfolder,[opt.prefix ff '.xml']),QAS,'write'); %struct('QAS',QAS,'QAM',QAM)
      end

      clear Yi Ym Yo Yos Ybc
      clear Ywm Ygm Ycsf Ybg

  end

  if nargout>2, varargout{3} = cat_qa_warnings; end
  if nargout>1, varargout{2} = QAR; end
  if nargout>0, varargout{1} = QAS; end 

end
%=======================================================================
function def=defaults
  % default parameter 
  def.verb       = 2;         % verbose level    [ 0=nothing | 1=points | 2*=results ]
  def.write_csv  = 2;         % final cms-file [ 0=dont write |1=write | 2=overwrite ]
  def.write_xml  = 1;         % images base xml-file
  def.sortQATm   = 1;         % sort QATm output
  def.orgval     = 0;         % original QAM results (no marks)
  def.avgfactor  = 2;         % 
  def.prefix     = 'cat_';    % intensity scaled  image
  def.mprefix    = 'm';       % prefix of the preprocessed image
  def.process    = 3;         % used image [ 0=T1 | 1=mT1 | 2=avg | 3=both ] 
  def.calc_MPC   = 0;
  def.calc_STC   = 0;
  def.calc_MJD   = 0;
  def.method     = 'spm';
  def.snspace    = [70,7,3];
  def.nogui      = exist('XT','var');
  def.MarkColor = cat_io_colormaps('marks+',40); 
end

function noise = estimateNoiseLevel(Ym,YM,r,rms,vx_vol)
% ----------------------------------------------------------------------
% noise estimation within Ym and YM.
% ----------------------------------------------------------------------
  if ~exist('vx_vol','var');
    vx_vol=[1 1 1]; 
  end
  if ~exist('r','var');
    r = 1;
  else
    r = min(10,max(max(vx_vol),r));
  end
  if ~exist('rms','var')
    rms = 1;
  end
  
  Ysd   = cat_vol_localstat(single(Ym),YM,r,4);
  noise = cat_stat_nanstat1d(Ysd(YM).^rms,'median').^(1/rms); 
end
%=======================================================================
