function cat_io_senderrormail(job)
% ______________________________________________________________________
% Function that prepares a mail for erros of CAT. 
% As input a set of directories is required. 
% It extract XML and PDF files of failed preprocessings.
% ______________________________________________________________________
% Robert Dahnke 
% $Revision: 1135 $  $Date: 2017-06-07 11:04:11 +0200 (Mi, 07 Jun 2017) $

%#ok<*TRYNC>

  if ~exist('job','var'), job = struct(); end
  
  if ~isfield(job,'errdirs')
    job.errdirs = cellstr(spm_select(Inf,'dir','Select processing or error directories'));
  else
    job.errdirs = cellstr(job.errdirs);
  end
  if isempty(job.errdirs), fprinft('Nothing to do. No input directory! \n'); return; end
  
  def.verb        = 1; % display progress
  def.optimize    = 1; % only unique errors
  def.onlyerrdirs = 1; % only error directories
  job = cat_io_checkinopt(job,def);
  
  % program version
  [nam,rev_cat] = cat_version;
  [nam,rev_spm] = spm('Ver');

  % initial message  
  emailSubject = sprintf('CAT %s error',rev_cat); 
  mainBody     = sprintf('Hi Christian,\\n\\nI got a serious problem with CAT %s in SPM %s with MATLAB %4.0f%s under %s.\\n',...
    rev_cat , rev_spm, round(cat_io_matlabversion/10) , char(97 + mod(cat_io_matlabversion,2)) , computer); 

  
  %% check errdirs
  %  -------------------------------------------------------------------
  if ~isempty(job.errdirs) && ~isempty(job.errdirs{1})
    if job.verb, fprintf('\nSearch for subdirectories .. '); end
    % remove twice entries
    job.errdirs = unique(job.errdirs);

    % remove old error dirs part 1
    for i=numel(job.errdirs):-1:1
      if ~isempty(strfind(job.errdirs{i},sprintf('cat%serror',rev_cat)))
         job.errdirs(i) = [];
      end
    end

    % find all subdirs
    job.errdirs = cat_vol_findfiles(job.errdirs,'*',struct('dirs',1));

    % remove old error dirs part 2
    for i=numel(job.errdirs):-1:1
      if ~isempty(strfind(job.errdirs{i},sprintf('cat%serror',rev_cat)))
         job.errdirs(i) = [];
      end
    end

    % remove twice subdirs
    for i=numel(job.errdirs):-1:1
      job.errdirs2 = job.errdirs; job.errdirs2{i} = ' '; 
      if any(~cellfun('isempty',(strfind(job.errdirs2,job.errdirs{i}))))
         job.errdirs(i) = [];
      end
    end
    % remove non error directories
    if job.onlyerrdirs 
      del = cellfun('isempty',strfind(job.errdirs,[filesep 'err'])); 
      job.errdirs(del) = [];
    end
    if isempty(job.errdirs) 
      fprintf('No errors directories "err" found.\n');
    end
  end
  if job.verb, fprintf('%d error directories found\n',numel(job.errdirs)); end 
  
  
  %% check for files
  %  -------------------------------------------------------------------
  if ~isempty(job.errdirs) && ~isempty(job.errdirs{1})
    % get the xml-files of the error dirs
    if job.verb, fprintf('\nSearch XML files .. '); end
    files = cat_vol_findfiles(job.errdirs,'cat*.xml');
    
    % remove old error dirs 
    for i=numel(files):-1:1
      if ~isempty(strfind(files{i},sprintf('cat%serror',rev_cat)))
         files(i) = [];
      end
    end
    if job.verb, fprintf('%d files found',numel(files)); end
    
    
    % load xml files
    if job.verb, fprintf('\nRead XML files .. '); end
    xml = cat_io_xml(files,struct(),'read',1)';
    
    % remove other version or files without error message
    if ~isfield(xml,'error'), files = []; xml = []; end
    for i=numel(files):-1:1
      if (xml(i).software.version_cat~=str2double(rev_cat)) || isempty(xml(i).error)
        files(i) = [];
        xml(i)   = []; 
      end
    end
    
    %% short error message 
    error = cell(size(xml));
    for i=1:numel(files)
      for ie=3:numel(xml(i).error)
        xml(i).error{ie} = strrep(xml(i).error{ie},' ',''); 
        errs = textscan(xml(i).error{ie},'%s','delimiter','-');
        xml(i).error{ie} = sprintf('%s:%s',errs{1}{2},errs{1}{1}); 
     end
      error{i}  = sprintf(sprintf('%%s%s',repmat(' | %s',1,numel(xml(i).error)-1)),xml(i).error{:});
    end
    %% remove twice error messages
    if job.optimize
      if job.verb, fprintf('\nRemove twice errors .. '); end; 
      [error,errornri,errornum] = unique(error); 
      for i=1:numel(error)
        error{i}  = sprintf('%4dx:%s',sum(errornum==i),error{i}); 
      end
      files = files(errornri);
      xml   = xml(errornri);
    end
    if job.verb, fprintf('%d files remain',numel(files)); end
    if ispc
      for i=numel(files):-1:1
        if any(strfind(files{i},' '))
          files(i) = [];
        end
      end
    end
    
    
    %% find pdfs
    pdfs = cell(''); pdf = files;
    for i=1:numel(files)
      [pp,ff] = spm_fileparts(files{i});
      pdf{i} = cat_vol_findfiles(pp,strrep([ff '.pdf'],'cat_','catreport_'),struct('maxdepth',1,'chararr',1)); 
      if ~isempty(pdf), pdfs = [pdfs;pdf{i}]; end %#ok<AGROW>
    end
    allfiles = [files;pdfs]; 
    
    
  
    %% remove tailing path names from the file directories
    if job.verb, fprintf('\nRead path names .. '); end
    for xmli=1:numel(xml)
      fn = fieldnames(xml(1).filedata); 
      for fni=1:numel(fn)
        for ai=1:numel(job.errdirs),
          if ischar(xml(xmli).filedata.(fn{fni}))
            xml(xmli).filedata.(fn{fni}) = strrep(xml(xmli).filedata.(fn{fni}),strrep(job.errdirs{ai},'/err',''),['..',filesep]); % pathname cat
          elseif iscell(xml(xmli).filedata.(fn{fni}))
            for xi=1:numel(xml(xmli).filedata.(fn{fni}))
              xml(xmli).filedata.(fn{fni}){xi} = strrep(xml(xmli).filedata.(fn{fni}){xi},strrep(job.errdirs{ai},'/err',''),['..',filesep]); % pathname cat
            end
          end
        end
      end
    end
    
    
    %% copy into packing directory
    packdir = fullfile(job.errdirs{1},sprintf('cat%serror',rev_cat)); 
    if exist(packdir,'dir'); 
      try 
        backdirfiles = cat_vol_findfiles(packdir,'*');
        for i=1:numel(backdirfiles)
          delete(backdirfiles{i});
        end
      end
      try
        rmdir(packdir,'s');
      end
    end
    if ~exist(packdir,'dir'); mkdir(packdir); end
    %%
    allfiless = spm_str_manip(allfiles,'c');
    for i=1:numel(allfiles)
      subdir = fullfile(packdir,allfiless{i}); 
      if ~exist(subdir,'dir'); mkdir(subdir); end
      if exist(allfiles{i},'file')
        [pp,ff,ee] = spm_fileparts(allfiles{i});
        try
          switch ee
            case '.xml' % update xml file 
              cat_io_xml(fullfile(subdir,[ff ee]),xml(i));
            otherwise
              if exist(allfiles{i},'file') 
                copyfile(allfiles{i},subdir,'f');
              end
          end
        catch
          fprintf('Failed to copy file "%s".\n',allfiles{i})
        end
      end
    end
    
    
    %% create ZIP file and add text
    if job.verb, fprintf('\nCreate ZIP file '); end
    zipfile = fullfile(job.errdirs{1},sprintf('cat%serror.zip',rev_cat)); 
    if exist(zipfile,'file'), delete(zipfile); end
    zip(zipfile,packdir);
    
    mainBody = sprintf('%sI send you the XML and PDF error files packed in the ZIP achive.\\n',mainBody);
    mainBody = sprintf('%s\\nBest regards\\n\\n',mainBody);
    
    mainBody = sprintf('%s\\n\\n    YOU HAVE TO ATTACH THE ZIP FILE BY YOURSELF. YOU CAN FIND IT IN: \\n\\n    %s\\n\\n',mainBody,strrep(zipfile,'\','\\'));
    
    mainBody = sprintf('%s\\n\\n\\nContent off the ZIP file:\\n',mainBody);
    
    filess = sort(spm_str_manip(files,'c'));
    for i=1:numel(files)
      for ai=1:numel(job.errdirs),
        filess{i} = strrep(strrep(filess{i},job.errdirs{ai},''),'\','\\');
      end
      if ~isempty(pdf{i}),  filess{i} = [filess{i} ' + pdf']; end
        
      mainBody = sprintf('%s  ..%s\\n',mainBody,filess{i});
    end
    
    %% add compressed error messages
    error = sort(error); 
    mainBody = sprintf('%s\\n\\nCompressed error messages:\\n%s',mainBody,sprintf('%s\\n',error{:}));
   
  else
    mainBody = sprintf('%s\\nBest regards\\n\\n\\n',mainBody);
  end
  
  % final creation of the mail 
  recipients = 'vbmweb@gmail.com';
  
  % 
  mainBodyHTML = strrep(mainBody,'\n','<br />');
  mainBodyMac  = sprintf(mainBody);
  mainBodyWinD = sprintf(mainBody);
  mainBodyWinD = strrep(mainBodyWinD,'\\\\','\');
  mainBodyWinD = strrep(mainBodyWinD,'\\','\');
  
  % windows
  mainBody = strrep(mainBody,' ','%20');
  mainBody = strrep(mainBody,'\n','%0A');
  mainBody = strrep(mainBody,'!','%21');
  mainBody = strrep(mainBody,'#','%23');
  mainBody = strrep(mainBody,'%%','%25');
  mainBody = strrep(mainBody,'*','%2A');
  mainBody = strrep(mainBody,'/','%2F');
  mainBody = strrep(mainBody,'<','%3C');
  mainBody = strrep(mainBody,'>','%3E');
  mainBody = strrep(mainBody,'?','%3F');
  mainBody = strrep(mainBody,'\\\\\\\\','\');
  mainBody = strrep(mainBody,'\\\\','\');
  mainBody = strrep(mainBody,'\\','\');
  
  
  %% print message, if mail does not work
 
  
  % create mail
  switch computer
    case {'MACI64','MACI32'}
      % here we need the " for each field
      fprintf('\nBEGIN MESSAGE:\n\n%s\n\nEND MESSAGE.\n\n',mainBodyMac);
      web(sprintf('mailto:"%s"?subject="%s"&body="%s"',recipients,emailSubject,mainBodyMac));
    case {'PCWIN','PCWIN64'}
      % it look strange, but it required only one " at the begining
      fprintf('\nBEGIN MESSAGE:\n\n%s\n\nEND MESSAGE.\n\n',mainBodyWinD);
      web(sprintf('mailto:"%s"?subject="%s"&body="%s"',recipients,emailSubject,mainBody),'-browser');
    otherwise %case {'GLNXA64'}
      %%
      fprintf('\nBEGIN MESSAGE:\n\n%s\n\nEND MESSAGE.\n\n',mainBodyMac);
      web(sprintf(['text://<head>' ...
        '<title>CAT12 - Error message</title>' ...
        '<meta http-equiv="Content-Type" content="text/html; charset=iso-8859-1" />' ...
        '</head>' ...
        '<body style="font: 10; color: #444;margin: 0px;padding: 0px;background: #F5F5F5;">' ...
        '<div id="header" style="color: white; background: black url(''%s'') right;' ...
      	'background-repeat: no-repeat;	width: 100%%; position: relative; height: 200px; ' ...
        'padding: 85px 0px 0px 20px; margin: 0px 0px 0px 0px; background-repeat: no-repeat">' ...
        '<br /><br /><h1>CAT Error Report</h1>' ...
        '</div>' ...
        '<div style="padding: 5px 20px 20px 20px;">' ... 
        '<p></p><p>The mailto action of the matlab WEB function does not work under LINUX. <br />' ...
        'Pleace copy and paste the default error message and the ZIP file by yourself</p>' ...
        '<h2>Recipient:</h2>%s</p>' ...
        '<h2>Email subject:</h2>%s</p>' ...
        '<h2>Message:</h2>' ...
        '<p>%s</p>' ...
        '</div></body>'],...
        fullfile(spm('dir'),'toolbox','cat12','html','images','contact.jpg'),...
        recipients,emailSubject,mainBodyHTML), '-new', '-notoolbar') ;
      
      
  end
  
  fprintf('Create CAT error report finished.\n')
end