function cat_io_xml2csv(PD,FN,fname)
% ______________________________________________________________________
%
% Function to convert a set of cat*.xml-files to a common csv-file.
%
%   cat_io_xml2csv([PD,FN])
%
%   PD = cell list of files or directories 
%   FN = cell list of fieldnames to convert (default = '*')
%
% Examples:
%   cat_io_xml2csv({'/my_cat_dir'},{'*'});
%
%   cat_io_xml2csv({'/my_cat_dir/control','/my_cat_dir/patients'},...
%                  {'qa.SM.vol_rel_CGW'});
%
% ______________________________________________________________________
% Robert Dahnke
% Structural Brain Mapping Group
% University Jena
%
% $Id: cat_io_xml2csv.m 767 2015-11-17 16:00:13Z gaser $
%
%#ok<*WNOFF,*WNON,*ASGLU>


  % get/check files
  if ~exist('PD','var') || isempty(PD) || strcmp(PD,'dirs')
    D = cellstr(spm_select([1 inf],'dir','Select directories with cat*.xml files',{},pwd));
    P = cell(size(D)); % list of cat-xml-files
  elseif strcmp(PD,'files')
    P = cellstr(spm_select([1 inf],'files','Select cat*.xml files',{},pwd,'cat_*.xml'));
    D = cell(size(P)); % list of directories
  else
    PD = cellstr(PD);
    P = cell(size(PD)); % list of cat-xml-files
    D = cell(size(PD)); % list of directories
    for di=1:numel(PD)
      if exist(PD{di},'dir')
        D{di} = PD{di};
      elseif exist(PD{di},'file')
        P{di}{1} = PD{di};
      end
    end
    clear PD;
  end
  % if directories are given, then find all cat-xml-files within
  if ~isempty(D)
    for di=numel(D):-1:1
      P{di}=cat_vol_findfiles(D{di},'cat_*.xml');
      if isempty(P{di})
        P(di)=[];
      end
    end
  end
  clear D;
  if isempty(P)
    fprintf('cat_io_xml2csv: No ''cat_*.xml''-files found!\n');
  else
    PP={};
    for di=numel(P)
      PP=[PP;P{di}]; %#ok<AGROW>
    end
    P=PP; clear PP;
  end
 
  
  % list of catDB-xml-files, if available
  Pdb = cell(size(P)); 
  for fi=1:numel(P)
    [pp,ff,ee] = spm_fileparts(P{fi});
    Pff{fi} = ff;
    Pdb{fi} = fullfile(pp,['catDB_' ff(5:end) ee]);
    if ~exist(Pdb{fi},'file')
      clear Pdb; 
      break
    end
  end  
  
  
  % load xml
  xml   = cat_io_xml(P);
  if exist('Pdb','var')
    xmldb = cat_io_xml(Pdb);
  end
  
  % fieldlist
  if ~exist('FN','var') || isempty(FN)
  %  FN = {'*'};
  %elseif ischar(FN) && strcmp(FN,'default')
    [FN,FNqm,FNdb] = cat_io_xml2csv_defaultsfields;
    if 0
      FN=[FN;FNqm];
    end
  end
  
  
  % filename
  if ~exist('fname','var') || isempty(fname)
    fname = fullfile(pwd,'cat12.csv');
  end
  
  
  % create table
  [xmlH,xmlT] = cat_io_struct2table(xml,FN);
  if exist('Pdb','var')
    [xmlHdb,xmlTdb] = cat_io_struct2table(xmldb,FNdb);
    xmlH = [xmlH, xmlHdb];
    xmlT = [xmlT, xmlTdb];
  end
  
  
  % save as csv-file
  %cat_io_csv(fname,['xml-file' xmlH; Pff' xmlT]);
  cat_io_csv(fname,[xmlH;xmlT]);

end
function [FNqa,FNqm,FNdb] = cat_io_xml2csv_defaultsfields
  FNdb = {
    'subject.Project';'subject.Group';'subject.Site';'subject.Subject';
    'subject.Age';'subject.Sex';'subject.Hand';'subject.DOB';'subject.DOS';
   %'subject.SD';'subject.RS';
    'subject.Weight';%'subject.Heigh';
   %'subject.MMSE';'subject.CDR';
   %'subject.Education;'subject.Medication;'subject.SES';
    'scanner.Manufacturer';'scanner.Fieldstr';'scanner.Model';
    'scanner.Protocol';'scanner.FlipAngle';'scanner.Plane';
    'scanner.TE';'scanner.TR';'scanner.TI'
    };
  FNqa = {
    'qa.FD.file';
    'qa.QMo.NCR';'qa.QMm.NCR';
    'qa.QMo.ICR';'qa.QMm.ICR';
    'qa.QMo.res_vx_vol';
    'qa.QMo.res_vol';
    'qa.QMo.res_isotropy';
    'qa.QMo.res_RMS';
    %'qa.QMo.NERR';'qa.QMm.NERR';
    %'qa.QMo.STC';'qa.QMm.STC';
    'qa.QMo.MPC';'qa.QMm.MPC';
    'qa.QMo.MJD';'qa.QMm.MJD';
    %'qa.QMo.contrast';'qa.QMm.contrast';
    %'qa.QMo.res_BB';  
    'qa.QMo.tissue_mn'; 
    'qa.QMo.tissue_std';
    'qa.SM.vol_TIV';'qa.SM.vol_rel_CGW';'qa.SM.vol_abs_CGW';
    }; 
  FNqm = setdiff(FNqa,{'qa.FD.file';});
  for i=1:numel(FNqm), FNqm{i}=strrep(FNqm{i},'qa','qam'); end 
end
