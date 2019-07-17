function varargout = cat_io_csv(filename,varargin)
% ______________________________________________________________________
% Writes and read csv-files (with sheet-like subparts) of a cell-array C 
% of chars and numbers. The action 'write' or 'read' is definend by the 
% output of this function - no output means write, whereas an output 
% will stand for read.
%
%   cat_io_csv(filename,C[,sheet,pos,opt])
%   C = cat_io_csv(filename[,sheet,pos,opt])
%
%   filename            = string with or without csv
%   C                   = cell with chars and numbers {'Hallo' 'Welt'; 1 2.3}
%   sheet               = '' - NOT WORKING YET
%   pos                 = 'A3:B4' - position of the Data
%   opt.delimiter       = ';'
%      .komma           = ',' 
%      .linedelimiter   = '\n' 
%      .format          = '%0.4f'
%
% Examples:
%   cat_io_csv('test',{'Hallo','Welt';1,2.4})
%   C=cat_io_csv('test.csv','','A1:C1')
%   cat_io_csv('test',rand(3,3),'','',struct('delimiter',',','komma','.'))
%
% TODO:
%   * adding of sheets to save different cells in one file
% ______________________________________________________________________
% Robert Dahnke 
% Structural Brain Mapping Group
% University Jena 
%
% $Id: cat_io_csv.m 1017 2016-09-23 05:55:23Z dahnke $
% ______________________________________________________________________

  if nargout>0, action='r'; else action='w'; end
  if ~exist('filename','var'),    
    filename = spm_select([0 1],'xml','Select *.csv files',{},pwd,'.*.csv');
    if isempty(filename)
      if nargout>0, varargout{1}=cell(); end
      return;
    end
  else
    [pp,ff,ee] = fileparts(filename);
    if ~strcmp(ee,'.csv'); ee='.csv'; end
    filename = fullfile(pp,[ff ee]); 
  end     
  if strcmp(action,'w')
    if nargin < 1+(nargout==0), C = {}; else C = varargin{1}; end
    if ~isa(C,'cell'); C = num2cell(C); end
  end
  if nargin < 2+(nargout==0), sheet = '';       else sheet = varargin{2-(nargout>0)}; end
  if nargin < 3+(nargout==0), pos   = '';       else pos   = varargin{3-(nargout>0)}; end
  if nargin < 4+(nargout==0), opt   = struct(); else opt   = varargin{4-(nargout>0)}; end  
  
  def.delimiter       = ',';
  def.komma           = '.'; 
  def.linedelimiter   = '\n'; 
  def.format          = '%0.4f';
  def.finaldelimiter  = 0;
  
  opt = cat_io_checkinopt(opt,def);
  if opt.komma==',' && opt.komma == opt.delimiter, opt.delimiter = ';'; end

  switch action
    case {'write','w'}
      writecsv(filename,C,sheet,pos,opt); 
    case {'read','r'}
      varargout{1} = readcsv(filename,sheet,pos,opt);
    otherwise
      error('MATLAB:CAT_IO_CSV:unkown_action','Unknown action ''%s''',action);
  end
end

function C=readcsv(filename,sheet,pos,opt)
% __________________________________________________________________________________________________
% read data as xls if matlab xlswrite works else will try to use the cvs-files.
% __________________________________________________________________________________________________
  
  % set filename and load if it exist
  if ~exist(filename,'file'), fprintf('File "%s" does not exist./n',filename); C={}; return; end

  % read file and convert from string to cell
  fid = fopen(filename);
  mv = version; mvi = strfind(mv,'R');
 % if str2double(mv(mvi+1:mvi+4)) < 2015 % old ... str2double(mv(mvi+1:mvi+4)) > 2013 &&
 %   C1  = textscan(fid,'%q','delimiter',opt.linedelimiter,'BufSize',2^24); C1=C1{1};
 % else % new
    C1  = textscan(fid,'%q','delimiter',opt.linedelimiter); C1=C1{1};
 % end
  fclose(fid);
  
  % na toll... textscan entfernt den zweiten ", aber lässt den ersten
  % stehen... die haben doch den arsch offen... 
  % The matlab textscan removes the second " of quoted csv-entries. 
  % Try to find the next delimiter and to replace the ". 
  for i=1:size(C1,1)
    Quote=strfind(C1{i},'"'); 
    for qi=numel(Quote):-1:1
      Delim=strfind(C1{i}(Quote(qi):end),opt.delimiter) + Quote(qi) - 1;
      if ~isempty(Delim) && strcmp(C1{i}(1:Delim(1)-1),'"')
        C1{i}=[C1{i}(1:Delim(1)-1) '"' C1{i}(Delim(1):end)];
      end
    end
  end

 

  % es gibt hier ein problem mit 
  % - sonderzeichen wie äöü
  % - anderen sonderzeichen 
  C1  = strrep(C1,'Ã¤','ä');
  C1  = strrep(C1,'Ã¼','ü');
  C1  = strrep(C1,'Ã¶','ö');
  

  for i=1:size(C1,1),
    try
      
    %  if numel(strfind(C1{i},';'))>0
      if isempty(C1{i})
        C2{i} = '';  %#ok<AGROW>
      else
%         if str2double(mv(mvi+1:mvi+4)) > 2013 && str2double(mv(mvi+1:mvi+4)) < 2015
%           C2{i}=textscan(C1{i},'%q','delimiter',opt.delimiter,'BufSize',2^24)'; C2{i}=C2{i}{1}';  %#ok<AGROW>
%         else
          C2{i}=textscan(C1{i},'%q','delimiter',opt.delimiter)'; C2{i}=C2{i}{1}';  %#ok<AGROW>
%         end
        %if size(C2{i},2)~=size(C2{1}), fprintf('Line %4d: %4d - %4d\n',i,size(C2{i},2),size(C2{1},2)); end
      end
    % fprintf('%4.0f-%4.0f\n',numel(strfind(C1{i},';')),numel(C2{i}));
    catch  %#ok<CTCH>
      fprintf('WARNING:cat_io_csv:readcsv: Can''t read line %d!\n',i); C2{i}=cell(1,numel(C2{1})); %#ok<AGROW>
    end
  end
  C3=cell(size(C2,2),max(cellfun('size',C2,1)));
  for i=1:size(C2,2)
    for j=1:size(C2{i},2)
      C3{i,j}=C2{i}{j};
    end
  end

  % set colum and row...
  if isempty(pos), C=C3; else C=readC(C3,pos); end

  % if a field could be interpreted as a number, then convert it to a float 
  % ??? if there is a comma otherwise to integer???
  for i=1:numel(C), if ~isnan(str2double(C{i})), id=strfind(C{i},','); C{i}(id)='.'; C{i} = str2double(C{i}); end; end

end
function writecsv(filename,C,sheet,pos,opt)
% __________________________________________________________________________________________________
% write data as xls if matlab xlswrite works else it export cvs-files.
% __________________________________________________________________________________________________
  % check if xlswrite could work
  
  % set colum and row
  if isempty(pos), else C=readC(C,pos); end

  for i=1:numel(C)
    if ~isempty(C{i})
      if isnumeric(C{i}), 
%           if isnan(C{i}) || isinf(C{i}), C{i}=sprintf('%0.0f',C{i});
%           else C{i}=sprintf(sprintf('%%0.0f,%%0%0.0f.0f',opt.digit),fix(C{i}),abs(C{i}*10^opt.digit-fix(C{i})*10^opt.digit));
%           end
      end
      if ischar(C{i}),  id=regexp(C{i},['[' opt.delimiter opt.linedelimiter ']']); C{i}(id)=[]; end
    end
  end

  % read old file if there is one an merge the cells where C isn't
  % specified the old value still exist
  if exist(filename,'file');
    fid = fopen(filename); 
    MO  = textscan(fid,'%s','delimiter',opt.linedelimiter); MO = MO{1};
    fclose(fid);

    for i=size(MO,1), tmp=textscan(MO{i},'%s','delimiter',opt.delimiter); CO{i,:} = tmp{1}; end %#ok<AGROW>

    % set size
    CC=cell(max(size(C),size(CO))); CC(1:size(C ,1),1:size(C ,2))=C;  C =CC;
    CC=cell(max(size(C),size(CO))); CC(1:size(CO,1),1:size(CO,2))=CO; CO=CC;
    clear CC;

    % merge
    for i=find(cellfun('isempty',C)==0); CO(i)=C(i); end
  end

  M=cell(size(C,1),1);
  for i=1:size(C,1)
    for j=1:size(C,2)
      if ~isstruct(C{i,j}) && ~iscell(C{i,j})
        if C{i,j}==round(C{i,j})
          M{i}=[M{i} num2str(C{i,j},'%d') opt.delimiter];
        else
          switch opt.komma
            case '.',   M{i}=[M{i} num2str(C{i,j},opt.format) opt.delimiter];
            otherwise,  M{i}=[M{i} strrep(num2str(C{i,j},opt.format),'.',opt.komma) opt.delimiter]; 
          end
        end
      else
        M{i}=[M{i} 'ERR:' class(C) opt.delimiter];
      end
    end
    if opt.finaldelimiter
      M{i}=[M{i} opt.linedelimiter];
    else
      M{i}=[M{i}(1:end-1) opt.linedelimiter];
    end
  end
  M=cell2mat(M');

  if ~exist(fileparts(filename),'dir'), mkdir(fileparts(filename));end
  
  f=fopen(filename,'w'); if f~=-1, fprintf(f,M); fclose(f); end;

end
function [Cpos,ijpos]=readC(C,pos)
  i=strfind(pos,':'); if ~isempty(i), pos(i)=[]; end                       % remove double points
  tmp=textscan(pos,'%[^1234567890]%d'); colum=tmp{1}; row=tmp{2};          % separate colum and row in pos-string
  if size(row,1)>2, row(3:end,1)=[]; colum(3:end,1)=[]; end                % remove to positions if there are to many
  ijpos(:,1)=sort(cell2mat(base27dec(colum)));                                  % convert to ij-position
  ijpos(:,2)=sort(double(row));
  CX=cell(max([size(C,1),ijpos(:,2)']),max([size(C,2),ijpos(:,1)'])); CX(1:size(C,1),1:size(C,2))=C;
  Cpos=CX(ijpos(1,end):ijpos(end,end),ijpos(1,1):ijpos(end,1));
end
function d = base27dec(s)
% copied from xlswrite.m
%--------------------------------------------------------------------------
%   BASE27DEC(S) returns the decimal of string S which represents a number in
%   base 27, expressed as 'A'..'Z', 'AA','AB'...'AZ', and so on. Note, there is
%   no zero so strictly we have hybrid base26, base27 number system.
%
%   Examples
%       base27dec('A') returns 1
%       base27dec('Z') returns 26
%       base27dec('IV') returns 256
%--------------------------------------------------------------------------
  d=cell(numel(s));
  if iscell(s)
    for i=1:numel(s), d{i} = base27dec(s{i}); end
  else
    if length(s) == 1
       d = s(1) -'A' + 1;
    else
        cumulative = 0;
        for i = 1:numel(s)-1
            cumulative = cumulative + 26.^i;
        end
        indexes_fliped = 1 + s - 'A';
        indexes = fliplr(indexes_fliped);
        indexes_in_cells = num2cell(indexes);
        d = cumulative + sub2ind(repmat(26, 1,numel(s)), indexes_in_cells{:});
    end
  end
end 