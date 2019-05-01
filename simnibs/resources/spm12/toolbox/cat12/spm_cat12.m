function spm_cat12(varargin)
% ______________________________________________________________________
% CAT12 Toolbox wrapper to start CAT with different user modes or 
% default files.  Changing the user mode requires restarting of CAT and
% SPM.  The expert user mode allows to control further parameters and  
% semi-evaluated functions, whereas the developer mode contain parameter
% for internal tests and unsafe functions.
% 
%   cat12(action)
%   
%   CAT user modes:
%     action = ['default','expert','developer'] 
%
%   CAT default files for other species (in development):
%     action = ['oldwoldmonkeys'|'greaterapes']
%
%   CAT start with own default files:
%     action = 'select' 
%     action = 'mypath/cat_defaults_mydefaults'
%
% ______________________________________________________________________
% Christian Gaser, Robert Dahnke
% $Id: spm_cat12.m 1250 2017-12-20 16:17:28Z gaser $


rev = '$Rev: 1250 $';
global deffile;
global cprintferror;  % temporary, because of JAVA errors in cat_io_cprintf ... 20160307
%try clearvars -global deffile;  end %#ok<TRYNC>

% start cat with different default file
catdir = fullfile(spm('dir'),'toolbox','cat12'); 
catdef = fullfile(catdir,'cat_defaults.m');
if nargin==0 && (isempty(deffile) || strcmp(deffile,catdef))
  deffile = catdef; 
  restartspm = 0;
elseif nargin==1 
  deffile = varargin{1}; 
  restartspm = 1;
elseif nargin==2
  deffile = varargin{1};
  catdef  = varargin{2};
  restartspm = 1;
else
  deffile = catdef; 
  restartspm = 1;
end


% choose filesspecies 
speciesdisp = '';
switch cat_get_defaults('extopts.species')
  case {'select','human','default','expert','developer'}
    % nothing to do
  otherwise
    % load default to remove animal settings
    try clearvars -global cat; end %#ok<TRYNC>
    [deffile_pp,deffile_ff] = fileparts(catdef);
    oldwkd = cd; 
    cd(deffile_pp);
    try clearvars -global cat; end %#ok<TRYNC>
    clear cat;
    eval(deffile_ff);
    eval('global cat;'); 
    cd(oldwkd);
end
switch lower(deffile) 
  case 'select'
    deffile = spm_select(1,'batch','Select CAT default file!','',catdir);
    if isempty(deffile) 
      return
    end
  case {'default','human'}
    mycat  = cat_get_defaults; 
    mycat.extopts.expertgui = 0;
    restartspm = 1;
    deffile = catdef; 
  case 'expert'
    mycat  = cat_get_defaults; 
    mycat.extopts.expertgui = 1;
    restartspm = 1;
    deffile = catdef; 
  case 'developer'
    mycat  = cat_get_defaults; 
    mycat.extopts.expertgui = 2;
    restartspm = 1;
    deffile = catdef; 
  case {'greaterapes','lesserapes','oldworldmonkeys','newworldmonkeys','mammals','chimpanzees','dogs',...
        'greaterape' ,'lesserape' ,'oldworldmonkey' ,'newworldmonkey', 'mammal', 'chimpanzee' ,'dog'}
    switch lower(deffile)
      case {'greaterapes','greaterape'},          species = 'ape_greater';     speciesdisp = ' (greater apes)';
      %case {'lesserapes','lesserape'},            species = 'ape_lesser';      speciesdisp = ' (lesser apes)';
      case {'oldworldmonkeys','oldworldmonkey'},  species = 'monkey_oldworld'; speciesdisp = ' (oldworld monkeys)';
      %case {'newworldmonkeys','newworldmonkey'},  species = 'monkey_newworld'; speciesdisp = ' (newworld monkeys)';
      %case {'mammals','mammal'},                  species = 'mammal';          speciesdisp = ' (mammal)';
      case {'chimpanzees','chimpanzee'},          species = 'chimpanzee';      speciesdisp = ' (chimpanzee)';
      case {'dogs','dog'},                        species = 'dog';             speciesdisp = ' (dogs)';
      otherwise
        error('CAT:unreadySpecies','Templates of species "%s" are not ready yet.\n',deffile);
    end  
    
    mycat                      = cat_get_defaults;
    % change TPM and user higher resolution and expect stronger bias
    mycat.opts.tpm             = {fullfile(spm('dir'),'toolbox','cat12','templates_animals',[species '_TPM.nii'])};
    mycat.opts.biasreg         = 0.001;                                  % less regularisation 
    mycat.opts.biasfwhm        = 50;                                     % stronger fields 
    mycat.opts.samp            = 2;                                      % smaller resampling
    % use species specific templates, higher resolution, stronger corrections and less affine registration (by SPM) 
    mycat.extopts.species      = species;  
    mycat.extopts.brainscale   = 200; % non-human brain volume in cm3 (from literature) or scaling in mm (check your data)
    mycat.extopts.darteltpm    = {fullfile(spm('dir'),'toolbox','cat12','templates_animals',[species '_Template_1.nii'])}; % Indicate first Dartel template
    mycat.extopts.shootingtpm  = {fullfile(spm('dir'),'toolbox','cat12','templates_animals',[species '_Template_0_GS.nii'])}; % Indicate first Shooting template
    mycat.extopts.cat12atlas   = {fullfile(spm('dir'),'toolbox','cat12','templates_animals',[species '_cat.nii'])};        % VBM atlas with major regions for VBM, SBM & ROIs
    mycat.extopts.brainmask    = {fullfile(spm('dir'),'toolbox','cat12','templates_animals',[species '_brainmask.nii'])};  % brainmask for affine registration
    mycat.extopts.T1           = {fullfile(spm('dir'),'toolbox','cat12','templates_animals',[species '_T1.nii'])};         % T1 for affine registration
    mycat.extopts.sanlm        = 2;                                     % ISARNLM for stronger corrections
    mycat.extopts.restype      = 'best';        
    mycat.extopts.resval       = [0.50 0.30];                           % higher interal resolution 
    mycat.extopts.APP          = 5;                                     % less affine registration, but full corrections (by SPM)
    mycat.extopts.vox          = 1.00;                                  % voxel size for normalized data 
    mycat.extopts.bb           = [[-inf -inf -inf];[inf inf inf]];      % template default
    mycat.extopts.expertgui    = 2;                                     % set to expert later ...
    mycat.extopts.ignoreErrors = 1;  
    switch species
      case 'monkey_oldworld'
        mycat.extopts.atlas = { ... 
          fullfile(spm('dir'),'toolbox','cat12','templates_animals','monkey_oldworld_atlas_inia19NeuroMaps.nii') 1 {'csf','gm','wm'} 1; 
          };
      case 'chimpanzee'
        mycat.extopts.atlas = { ... 
          fullfile(spm('dir'),'toolbox','cat12','templates_animals','chimpanzee_atlas_davi.nii') 1 {'csf','gm','wm'} 1; 
          };
      otherwise
        mycat.extopts.atlas = {}; 
        mycat.output.ROI    = 0;
    end
    
    restartspm = 1;
    deffile    = catdef; 
  otherwise
    % lazy input - no extension 
    [deffile_pp,deffile_ff,deffile_ee] = fileparts(deffile);
    if isempty(deffile_ee)
      deffile_ee = '.m';
    end
    % lazy input - no directory
    if isempty(deffile_pp) 
      if exist(fullfile(pwd,deffile_ff,deffile_ee),'file') 
        deffile_pp = pwd; 
      else
        deffile_pp = fullfile(spm('dir'),'toolbox','cat12'); 
      end
    end
    deffile = fullfile(deffile_pp,[deffile_ff,deffile_ee]); 

    if isempty(deffile) || ~exist(deffile,'file')
      help spm_cat12;
      error('CAT:unknownDefaultFile','Unknown action or nonexisting default file "%s".\n',deffile);
    end
end


% The cat12 global variable is created and localy destroyed, because we want to call the cat12 function. 
if exist('mycat','var') 
  try clearvars -global cat; end %#ok<TRYNC>
  eval('global cat; cat = mycat;'); 
else
  % set other defaultfile
  oldwkd = cd; 
  cd(deffile_pp);
  try clearvars -global cat; end %#ok<TRYNC>
  clear cat;
  eval(deffile_ff);
  eval('global cat;'); 
  cd(oldwkd);
end

% initialize SPM 
eval('global defaults;'); 
if restartspm 
  clear defaults; 
  spm_jobman('initcfg');
end
clear cat;

% temporary, because of JAVA errors in cat_io_cprintf ... 20160307
if cat_get_defaults('extopts.expertgui')<2
  cprintferror=1;
end


spm('FnBanner',mfilename,cat_version);
[Finter,Fgraph] = spm('FnUIsetup','CAT12.1');
url = fullfile(spm('Dir'),'toolbox','cat12','html','cat.html');
spm_help('!Disp',url,'',Fgraph,'Computational Anatomy Toolbox for SPM12');

[ST, RS] = cat_system('CAT_DumpCurv -h');
% because status will not give 0 for help output we have to check whether we can find the
% keyword "Usage" in output
if isempty(strfind(RS,'Usage'));
  if ispc
    [ST, RS] = system('systeminfo.exe');
  else
    [ST, RS] = system('uname -a');
  end
  cat_io_cmd(sprintf('\nWARNING: Surface processing will not work because\n(1) CAT binaries are not compatible to your system or\n(2) Antivirus software is blocking to execute binaries:\n%s\n',RS),'warning');
  fprintf('\n\nFor future support of your system please send this message to christian.gaser@uni-jena.de\n\n');
end


%% add some directories 
spm_select('PrevDirs',{fullfile(spm('dir'),'toolbox','cat12')});

%% command line output
cat_io_cprintf('silentreset');
switch cat_get_defaults('extopts.expertgui')
  case 0, expertguitext = '';
  case 1, expertguitext = ['Expert Mode' speciesdisp];
  case 2, expertguitext = ['Developer Mode' speciesdisp];
end
cat_io_cprintf([0.0 0.0 0.5],sprintf([ ...
    '\n' ...
    '   _______  ___  _______    \n' ...
    '  |  ____/ / _ \\\\ \\\\_   _/   ' expertguitext '\n' ...
    '  | |___  / /_\\\\ \\\\  | |     Computational Anatomy Toolbox\n' ...
    '  |____/ /_/   \\\\_\\\\ |_|     CAT12.1 - http://www.neuro.uni-jena.de\n\n']));
cat_io_cprintf([0.0 0.0 0.5],' CAT default file:\n\t%s\n\n',deffile); 

% call GUI
cat12('fig'); 
  

