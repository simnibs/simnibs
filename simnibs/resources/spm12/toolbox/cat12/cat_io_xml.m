function varargout = cat_io_xml(file,varargin)
% ______________________________________________________________________
% Import/export of a matlab structure from/to a xml file. Use functions 
% from Jaroslaw Tuszynski on MATLAB Central (xml_io_tools_2010_11_05). 
% Because the XML decoding is very slow, the data is further stored and 
% loaded (if available) as MATLAB MAT-file.
% 
%   cat_io_xml(file,S)      export structure to a xml-file
%   S = cat_io_xml(file)    import structure from a xml-file
%
% ______________________________________________________________________
% Copyright (c) 2007, Jaroslaw Tuszynski
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without 
% modification, are permitted provided that the following conditions are 
% met:
% 
%     * Redistributions of source code must retain the above copyright 
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright 
%       notice, this list of conditions and the following disclaimer in 
%       the documentation and/or other materials provided with the distribution
%       
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
% POSSIBILITY OF SUCH DAMAGE.
% ______________________________________________________________________
% $Id: cat_io_xml.m 1212 2017-11-10 12:39:37Z dahnke $


% Further comments:
% ______________________________________________________________________
% I also tried the "struct2xml" and "xml2struct" functions of Wouter 
% Falkena, but there were some problems for multi-element structures, 
% and fields where set as struct too and the only char is supported for
% datastoring.
% ______________________________________________________________________


  verbose = 0;
  if usejava('jvm')==0 || isdeployed
    warning('MATLAB:SPM:CAT:cat_io_xml:javaerror', ...
      'CAT-ERROR: CAT XML-im/export requires JVM! Read/Write only MAT file.\n');
    %varargout = {};
    %return;
  end
  if ~exist('file','var'),
    file = spm_select(Inf,'xml','Select *.xml files',{},pwd,'^cat.*.xml');
    if isempty(file)
      if nargout>0, varargout{1}=struct(); end
      return;
    end
  end
  if exist('varargin','var') 
    if numel(varargin)==0
      action='read';
    elseif numel(varargin)==1 
      if   ischar(varargin{1}), action='read'; % can only be read yet
      else S=varargin{1}; action='write';
      end
    elseif numel(varargin)==2
      S=varargin{1}; action=varargin{2};
      if ~isstruct(S)
        error('MATLAB:cat_io_xml','ERROR: Second input should be a structure.\n'); 
      end
    elseif numel(varargin)==3
      S=varargin{1}; action=varargin{2}; verbose=varargin{3};
      if ~isstruct(S)
        error('MATLAB:cat_io_xml','ERROR: Second input should be a structure.\n'); 
      end
    else  
      error('MATLAB:cat_io_xml','ERROR: To many inputs.\n');
    end
  end
  if (iscell(file) && numel(file)>100) || (ischar(file) && size(file,1)>100) 
    verbose = 1; 
  end
  
  % multi-file read 
  if strcmp(action,'read')
    varargout{1} = struct();
    if verbose, fprintf('% 6d/% 6d',0,numel(file)); end
    if iscell(file) && numel(file)>1 
      spm_progress_bar('Init',numel(file),...
        sprintf('read XML\n%d',numel(file)),'Files Completed'); 
      for fi=1:numel(file)
        try
          tmp = cat_io_xml(file{fi});
          fn = fieldnames(tmp);
          for fni = 1:numel(fn)
            varargout{1}(fi).(fn{fni}) = tmp.(fn{fni});
          end
          clear tmp;
        end
        if verbose, fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b% 6d/% 6d',fi,numel(file)); end
        spm_progress_bar('Set',fi);
      end

      spm_progress_bar('Clear');
      if verbose, fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b             \b\b\b\b\b\b\b\b\b\b\b\b\b'); end
      return
   
    elseif ischar(file) && size(file,1)>1
      spm_progress_bar('Init',size(file,1),...
        sprintf('read XML\n%s',size(file,1)),'Files Completed'); 
      for fi=1:numel(file)
        try
          tmp = cat_io_xml(file(fi,:));
          fn = fieldnames(tmp);
          for fni = 1:numel(fn)
            varargout{1}(fi).(fn{fni}) = tmp.(fn{fni});
          end
          clear tmp;
        end
        if verbose, fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b% 6d/% 6d',fi,numel(file)); end
        spm_progress_bar('Set',fi);
      end

      spm_progress_bar('Clear');
      if verbose, fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b             \b\b\b\b\b\b\b\b\b\b\b\b\b'); end
      return
    end
  end
  
  if iscell(file) && size(file,1)<=1, file = char(file); end
  
  [pp,ff,ee] = fileparts(file); if ~strcmp(ee,'.xml'), file = [file '.xml']; end
  if isempty(ff), return; end
  
  mfile = [file(1:end-4) '.mat']; 
  switch action
    case 'write'
    % ------------------------------------------------------------------  
      try
        S=orderfields(S); 
        if ~exist(fileparts(mfile),'dir'), try mkdir(fileparts(mfile)); end; end
        save(mfile,'S');
      catch %#ok<*NASGU> % can write xml file??
        error('MATLAB:cat_io_xml:writeErr','Can''t write MAT-file ''%s''!\n',mfile);
      end
      try
        if usejava('jvm') && ~isdeployed
          xml_write(file,S);
        end
      catch %#ok<*NASGU> % can write xml file??
        error('MATLAB:cat_io_xml:writeErr','Can''t write XML-file ''%s''!\n',file);
      end
      
      
    case 'write+'
    % ------------------------------------------------------------------  
    % WARNING: THIS ACTION NEED MUCH MORE WORK!!! 
    % ------------------------------------------------------------------  
      SN = orderfields(S); 
      
      if exist(mfile,'file')
        load(mfile,'S');
      elseif exist(file,'file') && usejava('jvm') && ~isdeployed
        try
          S = xml_read(file);
        catch 
          error('MATLAB:cat_io_xml:write+ReadErr','Can''t read XML-file ''%s'' for update!\n',file);
        end
      else
        S = struct(); 
      end
      
      if numel(S)>1 || numel(S)>1,
        error('MATLAB:cat_io_xml:write','Not implemented yet!\n');
      end
      
      S=cat_io_updateStruct(S,SN);
      
      try
        if ~exist(fileparts(mfile),'dir'), try mkdir(fileparts(mfile)); end; end
        save(mfile,'S');
      catch 
        error('MATLAB:cat_io_xml:writeErr','Can''t write MAT-file ''%s''!\n',mfile);
      end 
      try
        xml_write(file,S);
      catch 
        error('MATLAB:cat_io_xml:writeErr','Can''t write XML-file ''%s''!\n',file);
      end 
    
      
    case 'read'
    % ------------------------------------------------------------------
    % 
      if exist(mfile,'file')
        try 
          load(mfile,'S');
        catch
          error('MATLAB:cat_io_xml:readErr','Can''t read MAT-file ''%s'' for update!\n',mfile);
        end
      elseif exist(file,'file') && usejava('jvm') && ~isdeployed
        try 
          warning off
          S = xml_read(file);
          warning on
        catch 
          error('MATLAB:cat_io_xml:readErr','Can''t read XML-file ''%s'' for update!\n',file);
        end
        if ~exist(mfile,'file')
          try
            if ~exist(fileparts(mfile),'dir'), try, mkdir(fileparts(mfile)); end; end
            save(mfile,'S');
          catch
             error('MATLAB:cat_io_xml:writeErr','Can''t write MAT-file ''%s''!\n',mfile);
          end
        end
      elseif exist(file,'file') && ~(usejava('jvm') && ~isdeployed)
        S = struct(); 
      else
        error('MATLAB:cat_io_xml','"%s" does not exist!\n',file);
      end
    otherwise 
      error('MATLAB:cat_io_xml:read','Unknown action ''%s''!\n',action');
  end   
  
  if nargout==1, varargout{1} = S; end
          
end


function [tree, RootName, DOMnode] = xml_read(xmlfile, Pref)
%XML_READ reads xml files and converts them into Matlab's struct tree.
%
% DESCRIPTION
% tree = xml_read(xmlfile) reads 'xmlfile' into data structure 'tree'
%
% tree = xml_read(xmlfile, Pref) reads 'xmlfile' into data structure 'tree'
% according to your preferences
%
% [tree, RootName, DOMnode] = xml_read(xmlfile) get additional information
% about XML file
%
% INPUT:
%  xmlfile	URL or filename of xml file to read
%  Pref     Preferences:
%    Pref.ItemName - default 'item' - name of a special tag used to itemize
%                    cell arrays
%    Pref.ReadAttr - default true - allow reading attributes
%    Pref.ReadSpec - default true - allow reading special nodes
%    Pref.Str2Num  - default 'smart' - convert strings that look like numbers
%                   to numbers. Options: "always", "never", and "smart"
%    Pref.KeepNS   - default true - keep or strip namespace info
%    Pref.NoCells  - default true - force output to have no cell arrays
%    Pref.Debug    - default false - show mode specific error messages
%    Pref.NumLevels- default infinity - how many recursive levels are
%      allowed. Can be used to speed up the function by prunning the tree.
%    Pref.RootOnly - default true - output variable 'tree' corresponds to
%      xml file root element, otherwise it correspond to the whole file.
%    Pref.CellItem - default 'true' - leave 'item' nodes in cell notation.
% OUTPUT:
%  tree         tree of structs and/or cell arrays corresponding to xml file
%  RootName     XML tag name used for root (top level) node.
%               Optionally it can be a string cell array storing: Name of
%               root node, document "Processing Instructions" data and
%               document "comment" string
%  DOMnode      output of xmlread
%
% DETAILS:
% Function xml_read first calls MATLAB's xmlread function and than
% converts its output ('Document Object Model' tree of Java objects)
% to tree of MATLAB struct's. The output is in format of nested structs
% and cells. In the output data structure field names are based on
% XML tags, except in cases when tags produce illegal variable names.
%
% Several special xml node types result in special tags for fields of
% 'tree' nodes:
%  - node.CONTENT - stores data section of the node if other fields are
%    present. Usually data section is stored directly in 'node'.
%  - node.ATTRIBUTE.name - stores node's attribute called 'name'.
%  - node.COMMENT - stores node's comment section (string). For global
%    comments see "RootName" output variable.
%  - node.CDATA_SECTION - stores node's CDATA section (string).
%  - node.PROCESSING_INSTRUCTIONS - stores "processing instruction" child
%    node. For global "processing instructions" see "RootName" output variable.
%  - other special node types like: document fragment nodes, document type
%   nodes, entity nodes, notation nodes and processing instruction nodes
%   will be treated like regular nodes
%
% EXAMPLES:
%   MyTree=[];
%   MyTree.MyNumber = 13;
%   MyTree.MyString = 'Hello World';
%   xml_write('test.xml', MyTree);
%   [tree treeName] = xml_read ('test.xml');
%   disp(treeName)
%   gen_object_display()
%   % See also xml_examples.m
%
% See also:
%   xml_write, xmlread, xmlwrite
%
% Written by Jarek Tuszynski, SAIC, jaroslaw.w.tuszynski_at_saic.com
% References:
%  - Function inspired by Example 3 found in xmlread function.
%  - Output data structures inspired by xml_toolbox structures.

  % default preferences
  DPref.TableName  = {'tr','td'}; % name of a special tags used to itemize 2D cell arrays
  DPref.ItemName  = 'item'; % name of a special tag used to itemize 1D cell arrays
  DPref.CellItem  = false;  % leave 'item' nodes in cell notation
  DPref.ReadAttr  = true;   % allow reading attributes
  DPref.ReadSpec  = true;   % allow reading special nodes: comments, CData, etc.
  DPref.KeepNS    = true;   % Keep or strip namespace info
  DPref.Str2Num   = 'smart';% convert strings that look like numbers to numbers
  DPref.NoCells   = true;   % force output to have no cell arrays
  DPref.NumLevels = 1e10;   % number of recurence levels
  DPref.PreserveSpace = false; % Preserve or delete spaces at the beggining and the end of stings?
  RootOnly        = true;   % return root node  with no top level special nodes
  Debug           = false;  % show specific errors (true) or general (false)?
  tree            = [];
  RootName        = [];

  % Check Matlab Version
  v = ver('MATLAB');
  version = str2double(regexp(v.Version, '\d.\d','match','once'));
  if (version<7.1)
    error('Your MATLAB version is too old. You need version 7.1 or newer.');
  end

  % read user preferences
  if (nargin>1)
    if (isfield(Pref, 'TableName')), DPref.TableName = Pref.TableName; end
    if (isfield(Pref, 'ItemName' )), DPref.ItemName  = Pref.ItemName;  end
    if (isfield(Pref, 'CellItem' )), DPref.CellItem  = Pref.CellItem;  end
    if (isfield(Pref, 'Str2Num'  )), DPref.Str2Num   = Pref.Str2Num ;  end
    if (isfield(Pref, 'NoCells'  )), DPref.NoCells   = Pref.NoCells ;  end
    if (isfield(Pref, 'NumLevels')), DPref.NumLevels = Pref.NumLevels; end
    if (isfield(Pref, 'ReadAttr' )), DPref.ReadAttr  = Pref.ReadAttr;  end
    if (isfield(Pref, 'ReadSpec' )), DPref.ReadSpec  = Pref.ReadSpec;  end
    if (isfield(Pref, 'KeepNS'   )), DPref.KeepNS    = Pref.KeepNS;    end
    if (isfield(Pref, 'RootOnly' )), RootOnly        = Pref.RootOnly;  end
    if (isfield(Pref, 'Debug'    )), Debug           = Pref.Debug   ;  end
    if (isfield(Pref, 'PreserveSpace')), DPref.PreserveSpace = Pref.PreserveSpace; end
  end
  if ischar(DPref.Str2Num), % convert from character description to numbers
    DPref.Str2Num = find(strcmpi(DPref.Str2Num, {'never', 'smart', 'always'}))-1;
    if isempty(DPref.Str2Num), DPref.Str2Num=1; end % 1-smart by default
  end

  % read xml file using Matlab function
  if isa(xmlfile, 'org.apache.xerces.dom.DeferredDocumentImpl');
    % if xmlfile is a DOMnode than skip the call to xmlread
    try
      %try
        DOMnode = xmlfile;
      %catch ME
      %  error('Invalid DOM node: \n%s.', getReport(ME));
      %end
    catch %#ok<CTCH> catch for mablab versions prior to 7.5
      error('Invalid DOM node. \n');
    end
  else         % we assume xmlfile is a filename
    if (Debug) % in debuging mode crashes are allowed
      DOMnode = xmlread(xmlfile);
    else       % in normal mode crashes are not allowed
      try
        %try
          DOMnode = xmlread(xmlfile);
        %catch ME
        %  error('Failed to read XML file %s: \n%s',xmlfile, getReport(ME));
        %end
      catch %#ok<CTCH> catch for mablab versions prior to 7.5
        error('Failed to read XML file %s\n',xmlfile);
      end
    end
  end
  Node = DOMnode.getFirstChild;

  % Find the Root node. Also store data from Global Comment and Processing
  %  Instruction nodes, if any.
  GlobalTextNodes = cell(1,3);
  GlobalProcInst  = [];
  GlobalComment   = [];
  GlobalDocType   = [];
  while (~isempty(Node))
    if (Node.getNodeType==Node.ELEMENT_NODE)
      RootNode=Node;
    elseif (Node.getNodeType==Node.PROCESSING_INSTRUCTION_NODE)
      data   = strtrim(char(Node.getData));
      target = strtrim(char(Node.getTarget));
      GlobalProcInst = [target, ' ', data];
      GlobalTextNodes{2} = GlobalProcInst;
    elseif (Node.getNodeType==Node.COMMENT_NODE)
      GlobalComment = strtrim(char(Node.getData));
      GlobalTextNodes{3} = GlobalComment;
      %   elseif (Node.getNodeType==Node.DOCUMENT_TYPE_NODE)
      %     GlobalTextNodes{4} = GlobalDocType;
    end
    Node = Node.getNextSibling;
  end

  % parse xml file through calls to recursive DOMnode2struct function
  if (Debug)   % in debuging mode crashes are allowed
    [tree RootName] = DOMnode2struct(RootNode, DPref, 1);
  else         % in normal mode crashes are not allowed
    try
      %try
        [tree RootName] = DOMnode2struct(RootNode, DPref, 1);
      %catch ME
      %  error('Unable to parse XML file %s: \n %s.',xmlfile, getReport(ME));
      %end
    catch %#ok<CTCH> catch for mablab versions prior to 7.5
      error('Unable to parse XML file %s.',xmlfile);
    end
  end

  % If there were any Global Text nodes than return them
  if (~RootOnly)
    if (~isempty(GlobalProcInst) && DPref.ReadSpec)
      t.PROCESSING_INSTRUCTION = GlobalProcInst;
    end
    if (~isempty(GlobalComment) && DPref.ReadSpec)
      t.COMMENT = GlobalComment;
    end
    if (~isempty(GlobalDocType) && DPref.ReadSpec)
      t.DOCUMENT_TYPE = GlobalDocType;
    end
    t.(RootName) = tree;
    tree=t;
  end
  if (~isempty(GlobalTextNodes))
    GlobalTextNodes{1} = RootName;
    RootName = GlobalTextNodes;
  end
end

% -- begin of mathworks code that is equaly slow ... delete this later --
          function theStruct = parseXML(filename)
            % PARSEXML Convert XML file to a MATLAB structure.
            try
               tree = xmlread(filename);
            catch
               error('Failed to read XML file %s.',filename);
            end

            % Recurse over child nodes. This could run into problems 
            % with very deeply nested trees.
            try
               theStruct = parseChildNodes(tree);
            catch
               error('Unable to parse XML file %s.',filename);
            end
          end

          % ----- Local function PARSECHILDNODES -----
          function children = parseChildNodes(theNode)
            % Recurse over node children.
            children = [];
            if theNode.hasChildNodes
               childNodes = theNode.getChildNodes;
               numChildNodes = childNodes.getLength;
               allocCell = cell(1, numChildNodes);

               children = struct(             ...
                  'Name', allocCell, 'Attributes', allocCell,    ...
                  'Data', allocCell, 'Children', allocCell);

                for count = 1:numChildNodes
                    theChild = childNodes.item(count-1);
                    children(count) = makeStructFromNode(theChild);
                end
            end
          end

          % ----- Local function MAKESTRUCTFROMNODE -----
          function nodeStruct = makeStructFromNode(theNode)
            % Create structure of node info.

            nodeStruct = struct(                        ...
               'Name', char(theNode.getNodeName),       ...
               'Attributes', parseAttributes(theNode),  ...
               'Data', '',                              ...
               'Children', parseChildNodes(theNode));

            if any(strcmp(methods(theNode), 'getData'))
               nodeStruct.Data = char(theNode.getData); 
            else
               nodeStruct.Data = '';
            end
          end

          % ----- Local function PARSEATTRIBUTES -----
          function attributes = parseAttributes(theNode)
            % Create attributes structure.

            attributes = [];
            if theNode.hasAttributes
               theAttributes = theNode.getAttributes;
               numAttributes = theAttributes.getLength;
               allocCell = cell(1, numAttributes);
               attributes = struct('Name', allocCell, 'Value', ...
                                   allocCell);

               for count = 1:numAttributes
                  attrib = theAttributes.item(count-1);
                  attributes(count).Name = char(attrib.getName);
                  attributes(count).Value = char(attrib.getValue);
               end
            end
          end
          % -- end of mathworks code --


  % =======================================================================
  %  === DOMnode2struct Function ===========================================
  %  =======================================================================
  function [s TagName LeafNode] = DOMnode2struct(node, Pref, level)

    % === Step 1: Get node name and check if it is a leaf node ==============
    [TagName LeafNode] = NodeName(node, Pref.KeepNS);
    s = []; % initialize output structure

    % === Step 2: Process Leaf Nodes (nodes with no children) ===============
    if (LeafNode)
      if (LeafNode>1 && ~Pref.ReadSpec), LeafNode=-1; end % tags only so ignore special nodes
      if (LeafNode>0) % supported leaf node types
        try
          %try         % use try-catch: errors here are often due to VERY large fields (like images) that overflow java memory
            s = char(node.getData);
            if (isempty(s)), s = ' '; end                              % make it a string
            % for some reason current xmlread 'creates' a lot of empty text
            % fields with first chatacter=10 - those will be deleted.
            if (~Pref.PreserveSpace || s(1)==10) 
              if (isspace(s(1)) || isspace(s(end))), s = strtrim(s); end % trim speces is any
            end
            if (LeafNode==1), s=str2var(s, Pref.Str2Num, 0); end       % convert to number(s) if needed
          %catch ME    % catch for mablab versions 7.5 and higher
          %  warning('xml_io_tools:read:LeafRead', ...
          %    'This leaf node could not be read and was ignored. ');
          %  getReport(ME)
          %end
        catch         %#ok<CTCH> catch for mablab versions prior to 7.5
          warning('xml_io_tools:read:LeafRead', ...
            'This leaf node could not be read and was ignored. ');
        end
      end
      if (LeafNode==3) % ProcessingInstructions need special treatment
        target = strtrim(char(node.getTarget));
        s = [target, ' ', s];
      end
      return % We are done the rest of the function deals with nodes with children
    end
    if (level>Pref.NumLevels+1), return; end % if Pref.NumLevels is reached than we are done

    % === Step 3: Process nodes with children ===============================
    if (node.hasChildNodes)        % children present
      Child  = node.getChildNodes; % create array of children nodes
      nChild = Child.getLength;    % number of children

      % --- pass 1: how many children with each name -----------------------
      f = [];
      for iChild = 1:nChild        % read in each child
        [cname cLeaf] = NodeName(Child.item(iChild-1), Pref.KeepNS);
        if (cLeaf<0), continue; end % unsupported leaf node types
        if (~isfield(f,cname)),
          f.(cname)=0;           % initialize first time I see this name
        end
        f.(cname) = f.(cname)+1; % add to the counter
      end                        % end for iChild
      % text_nodes become CONTENT & for some reason current xmlread 'creates' a
      % lot of empty text fields so f.CONTENT value should not be trusted
      if (isfield(f,'CONTENT') && f.CONTENT>2), f.CONTENT=2; end

      % --- pass 2: store all the children as struct of cell arrays ----------
      for iChild = 1:nChild        % read in each child
        [c cname cLeaf] = DOMnode2struct(Child.item(iChild-1), Pref, level+1);
        if (cLeaf && isempty(c))   % if empty leaf node than skip
          continue;                % usually empty text node or one of unhandled node types
        elseif (nChild==1 && cLeaf==1)
          s=c;                     % shortcut for a common case
        else                       % if normal node
          if (level>Pref.NumLevels), continue; end
          n = f.(cname);           % how many of them in the array so far?
          if (~isfield(s,cname))   % encountered this name for the first time
            if (n==1)              % if there will be only one of them ...
              s.(cname) = c;       % than save it in format it came in
            else                   % if there will be many of them ...
              s.(cname) = cell(1,n);
              s.(cname){1} = c;    % than save as cell array
            end
            f.(cname) = 1;         % initialize the counter
          else                     % already have seen this name
            s.(cname){n+1} = c;    % add to the array
            f.(cname) = n+1;       % add to the array counter
          end
        end
      end   % for iChild
    end % end if (node.hasChildNodes)

    % === Step 4: Post-process struct's created for nodes with children =====
    if (isstruct(s))
      fields = fieldnames(s);
      nField = length(fields);

      % Detect structure that looks like Html table and store it in cell Matrix
      if (nField==1 && strcmpi(fields{1},Pref.TableName{1}))
        tr = s.(Pref.TableName{1});
        fields2 = fieldnames(tr{1});
        if (length(fields2)==1 && strcmpi(fields2{1},Pref.TableName{2}))
          % This seems to be a special structure such that for 
          % Pref.TableName = {'tr','td'} 's' corresponds to 
          %    <tr> <td>M11</td> <td>M12</td> </tr>
          %    <tr> <td>M12</td> <td>M22</td> </tr>
          % Recognize it as encoding for 2D struct
          nr = length(tr);
          for r = 1:nr
            row = tr{r}.(Pref.TableName{2});
            Table(r,1:length(row)) = row; %#ok<AGROW>
          end
          s = Table;
        end
      end

      % --- Post-processing: convert 'struct of cell-arrays' to 'array of structs'
      % Example: let say s has 3 fields s.a, s.b & s.c  and each field is an
      % cell-array with more than one cell-element and all 3 have the same length.
      % Then change it to array of structs, each with single cell.
      % This way element s.a{1} will be now accessed through s(1).a
      vec = zeros(size(fields));
      for i=1:nField, vec(i) = f.(fields{i}); end
      if (numel(vec)>1 && vec(1)>1 && var(vec)==0)  % convert from struct of
        s = cell2struct(struct2cell(s), fields, 1); % arrays to array of struct
      end % if anyone knows better way to do above conversion please let me know.

    end

    % === Step 5: Process nodes with attributes =============================
    if (node.hasAttributes && Pref.ReadAttr)
      if (~isstruct(s)),              % make into struct if is not already
        ss.CONTENT=s;
        s=ss;
      end
      Attr  = node.getAttributes;     % list of all attributes
      for iAttr = 1:Attr.getLength    % for each attribute
        name  = char(Attr.item(iAttr-1).getName);  % attribute name
        name  = str2varName(name, Pref.KeepNS);    % fix name if needed
        value = char(Attr.item(iAttr-1).getValue); % attribute value
        value = str2var(value, Pref.Str2Num, 1);   % convert to number if possible
        s.ATTRIBUTE.(name) = value;   % save again
      end                             % end iAttr loop
    end % done with attributes
    if (~isstruct(s)), return; end %The rest of the code deals with struct's

    % === Post-processing: fields of "s"
    % convert  'cell-array of structs' to 'arrays of structs'
    fields = fieldnames(s);     % get field names
    nField = length(fields);
    for iItem=1:length(s)       % for each struct in the array - usually one
      for iField=1:length(fields)
        field = fields{iField}; % get field name
        % if this is an 'item' field and user want to leave those as cells
        % than skip this one
        if (strcmpi(field, Pref.ItemName) && Pref.CellItem), continue; end
        x = s(iItem).(field);
        if (iscell(x) && all(cellfun(@isstruct,x(:))) && numel(x)>1) % it's cell-array of structs
          % numel(x)>1 check is to keep 1 cell-arrays created when Pref.CellItem=1
          try                           % this operation fails sometimes
            % example: change s(1).a{1}.b='jack'; s(1).a{2}.b='john'; to
            % more convinient s(1).a(1).b='jack'; s(1).a(2).b='john';
            s(iItem).(field) = [x{:}]';  %#ok<AGROW> % converted to arrays of structs
          catch %#ok<CTCH>
            % above operation will fail if s(1).a{1} and s(1).a{2} have
            % different fields. If desired, function forceCell2Struct can force
            % them to the same field structure by adding empty fields.
            if (Pref.NoCells)
              s(iItem).(field) = forceCell2Struct(x); %#ok<AGROW>
            end
          end % end catch
        end
      end
    end

    % === Step 4: Post-process struct's created for nodes with children =====

    % --- Post-processing: remove special 'item' tags ---------------------
    % many xml writes (including xml_write) use a special keyword to mark
    % arrays of nodes (see xml_write for examples). The code below converts
    % s.item to s.CONTENT
    ItemContent = false;
    if (isfield(s,Pref.ItemName))
      s.CONTENT = s.(Pref.ItemName);
      s = rmfield(s,Pref.ItemName);
      ItemContent = Pref.CellItem; % if CellItem than keep s.CONTENT as cells
    end

    % --- Post-processing: clean up CONTENT tags ---------------------
    % if s.CONTENT is a cell-array with empty elements at the end than trim
    % the length of this cell-array. Also if s.CONTENT is the only field than
    % remove .CONTENT part and store it as s.
    if (isfield(s,'CONTENT'))
      if (iscell(s.CONTENT) && isvector(s.CONTENT))
        x = s.CONTENT;
        for i=numel(x):-1:1, if ~isempty(x{i}), break; end; end
        if (i==1 && ~ItemContent)
          s.CONTENT = x{1};   % delete cell structure
        else
          s.CONTENT = x(1:i); % delete empty cells
        end
      end
      if (nField==1)
        if (ItemContent)
          ss = s.CONTENT;       % only child: remove a level but ensure output is a cell-array
          s=[]; s{1}=ss;
        else
          s = s.CONTENT;        % only child: remove a level
        end
      end
    end
  end
  % =======================================================================
  %  === forceCell2Struct Function =========================================
  %  =======================================================================
  function s = forceCell2Struct(x)
    % Convert cell-array of structs, where not all of structs have the same
    % fields, to a single array of structs

    % Convert 1D cell array of structs to 2D cell array, where each row
    % represents item in original array and each column corresponds to a unique
    % field name. Array "AllFields" store fieldnames for each column
    AllFields = fieldnames(x{1});     % get field names of the first struct
    CellMat = cell(length(x), length(AllFields));
    for iItem=1:length(x)
      fields = fieldnames(x{iItem});  % get field names of the next struct
      for iField=1:length(fields)     % inspect all fieldnames and find those
        field = fields{iField};       % get field name
        col = find(strcmp(field,AllFields),1);
        if isempty(col)               % no column for such fieldname yet
          AllFields = [AllFields; field]; %#ok<AGROW>
          col = length(AllFields);    % create a new column for it
        end
        CellMat{iItem,col} = x{iItem}.(field); % store rearanged data
      end
    end
    % Convert 2D cell array to array of structs
    s = cell2struct(CellMat, AllFields, 2);
  end
  % =======================================================================
  %  === str2var Function ==================================================
  %  =======================================================================
  function val=str2var(str, option, attribute)
    % Can this string 'str' be converted to a number? if so than do it.
    val = str;
    len = numel(str);
    if (len==0    || option==0), return; end % Str2Num="never" of empty string -> do not do enything
    if (len>10000 && option==1), return; end % Str2Num="smart" and string is very long -> probably base64 encoded binary
    digits = '(Inf)|(NaN)|(pi)|[\t\n\d\+\-\*\.ei EI\[\]\;\,]';
    s = regexprep(str, digits, ''); % remove all the digits and other allowed characters
    if (~all(~isempty(s)))          % if nothing left than this is probably a number
      if (~isempty(strfind(str, ' '))), option=2; end %if str has white-spaces assume by default that it is not a date string
      if (~isempty(strfind(str, '['))), option=2; end % same with brackets
      str(strfind(str, '\n')) = ';';% parse data tables into 2D arrays, if any
      if (option==1)                % the 'smart' option
        try                         % try to convert to a date, like 2007-12-05
          datenum(str);             % if successful than leave it as string
        catch                       %#ok<CTCH> % if this is not a date than ...
          option=2;                 % ... try converting to a number
        end
      end
      if (option==2)
        if (attribute)
          num = str2double(str);      % try converting to a single number using sscanf function
          if isnan(num), return; end  % So, it wasn't really a number after all    
        else
          num = str2num(str);         %#ok<ST2NM> % try converting to a single number or array using eval function
        end
        if(isnumeric(num) && numel(num)>0), val=num; end % if convertion to a single was succesful than save
      end
    elseif ((str(1)=='[' && str(end)==']') || (str(1)=='{' && str(end)=='}')) % this looks like a (cell) array encoded as a string
      try 
        val = eval(str); 
      catch              %#ok<CTCH>
        val = str; 
      end                     
    elseif (~attribute)   % see if it is a boolean array with no [] brackets
      str1 = lower(str);
      str1 = strrep(str1, 'false', '0');
      str1 = strrep(str1, 'true' , '1');
      s = regexprep(str1, '[01 \;\,]', ''); % remove all 0/1, spaces, commas and semicolons 
      if (~all(~isempty(s)))          % if nothing left than this is probably a boolean array
        num  = str2num(str1); %#ok<ST2NM>
        if(isnumeric(num) && numel(num)>0), val = (num>0);  end % if convertion was succesful than save as logical
      end
    end
  end
  % =======================================================================
  %  === str2varName Function ==============================================
  %  =======================================================================
  function str = str2varName(str, KeepNS)
    % convert a sting to a valid matlab variable name
    if(KeepNS)
      str = regexprep(str,':','_COLON_', 'once', 'ignorecase');
    else
      k = strfind(str,':');
      if (~isempty(k))
        str = str(k+1:end);
      end
    end
    str = regexprep(str,'-','_DASH_'  ,'once', 'ignorecase');
    if (~isvarname(str)) && (~iskeyword(str))
      str = genvarname(str);
    end
  end
  % =======================================================================
  %  === NodeName Function =================================================
  %  =======================================================================
  function [Name LeafNode] = NodeName(node, KeepNS)
    % get node name and make sure it is a valid variable name in Matlab.
    % also get node type:
    %   LeafNode=0 - normal element node,
    %   LeafNode=1 - text node
    %   LeafNode=2 - supported non-text leaf node,
    %   LeafNode=3 - supported processing instructions leaf node,
    %   LeafNode=-1 - unsupported non-text leaf node
    switch (node.getNodeType)
      case node.ELEMENT_NODE
        Name = char(node.getNodeName);% capture name of the node
        Name = str2varName(Name, KeepNS);     % if Name is not a good variable name - fix it
        LeafNode = 0;
      case node.TEXT_NODE
        Name = 'CONTENT';
        LeafNode = 1;
      case node.COMMENT_NODE
        Name = 'COMMENT';
        LeafNode = 2;
      case node.CDATA_SECTION_NODE
        Name = 'CDATA_SECTION';
        LeafNode = 2;
      case node.DOCUMENT_TYPE_NODE
        Name = 'DOCUMENT_TYPE';
        LeafNode = 2;
      case node.PROCESSING_INSTRUCTION_NODE
        Name = 'PROCESSING_INSTRUCTION';
        LeafNode = 3;
      otherwise
        NodeType = {'ELEMENT','ATTRIBUTE','TEXT','CDATA_SECTION', ...
          'ENTITY_REFERENCE', 'ENTITY', 'PROCESSING_INSTRUCTION', 'COMMENT',...
          'DOCUMENT', 'DOCUMENT_TYPE', 'DOCUMENT_FRAGMENT', 'NOTATION'};
        Name = char(node.getNodeName);% capture name of the node
        warning('xml_io_tools:read:unkNode', ...
          'Unknown node type encountered: %s_NODE (%s)', NodeType{node.getNodeType}, Name);
        LeafNode = -1;
    end
  end

  
function DOMnode = xml_write(filename, tree, RootName, Pref)
%XML_WRITE  Writes Matlab data structures to XML file
%
% DESCRIPTION
% xml_write( filename, tree) Converts Matlab data structure 'tree' containing
% cells, structs, numbers and strings to Document Object Model (DOM) node
% tree, then saves it to XML file 'filename' using Matlab's xmlwrite
% function. Optionally one can also use alternative version of xmlwrite
% function which directly calls JAVA functions for XML writing without
% MATLAB middleware. This function is provided as a patch to existing
% bugs in xmlwrite (in R2006b).
%
% xml_write(filename, tree, RootName, Pref) allows you to specify
% additional preferences about file format
%
% DOMnode = xml_write([], tree) same as above except that DOM node is
% not saved to the file but returned.
%
% INPUT
%   filename     file name
%   tree         Matlab structure tree to store in xml file.
%   RootName     String with XML tag name used for root (top level) node
%                Optionally it can be a string cell array storing: Name of
%                root node, document "Processing Instructions" data and
%                document "comment" string
%   Pref         Other preferences:
%     Pref.ItemName - default 'item' -  name of a special tag used to
%                     itemize cell or struct arrays
%     Pref.XmlEngine - let you choose the XML engine. Currently default is
%       'Xerces', which is using directly the apache xerces java file.
%       Other option is 'Matlab' which uses MATLAB's xmlwrite and its
%       XMLUtils java file. Both options create identical results except in
%       case of CDATA sections where xmlwrite fails.
%     Pref.CellItem - default 'true' - allow cell arrays to use 'item'
%       notation. See below.
%    Pref.RootOnly - default true - output variable 'tree' corresponds to
%       xml file root element, otherwise it correspond to the whole file.
%     Pref.StructItem - default 'true' - allow arrays of structs to use
%       'item' notation. For example "Pref.StructItem = true" gives:
%         <a>
%           <b>
%             <item> ... <\item>
%             <item> ... <\item>
%           <\b>
%         <\a>
%       while "Pref.StructItem = false" gives:
%         <a>
%           <b> ... <\b>
%           <b> ... <\b>
%         <\a>
%
%
% Several special xml node types can be created if special tags are used
% for field names of 'tree' nodes:
%  - node.CONTENT - stores data section of the node if other fields
%    (usually ATTRIBUTE are present. Usually data section is stored
%    directly in 'node'.
%  - node.ATTRIBUTE.name - stores node's attribute called 'name'.
%  - node.COMMENT - create comment child node from the string. For global
%    comments see "RootName" input variable.
%  - node.PROCESSING_INSTRUCTIONS - create "processing instruction" child
%    node from the string. For global "processing instructions" see
%    "RootName" input variable.
%  - node.CDATA_SECTION - stores node's CDATA section (string). Only works
%    if Pref.XmlEngine='Xerces'. For more info, see comments of F_xmlwrite.
%  - other special node types like: document fragment nodes, document type
%    nodes, entity nodes and notation nodes are not being handled by
%    'xml_write' at the moment.
%
% OUTPUT
%   DOMnode      Document Object Model (DOM) node tree in the format
%                required as input to xmlwrite. (optional)
%
% EXAMPLES:
%   MyTree=[];
%   MyTree.MyNumber = 13;
%   MyTree.MyString = 'Hello World';
%   xml_write('test.xml', MyTree);
%   type('test.xml')
%   %See also xml_tutorial.m
%
% See also
%   xml_read, xmlread, xmlwrite
%
% Written by Jarek Tuszynski, SAIC, jaroslaw.w.tuszynski_at_saic.com

  % Check Matlab Version
  v = ver('MATLAB');
  v = str2double(regexp(v.Version, '\d.\d','match','once'));
  if (v<7)
    error('Your MATLAB version is too old. You need version 7.0 or newer.');
  end

  % default preferences
  DPref.TableName  = {'tr','td'}; % name of a special tags used to itemize 2D cell arrays
  DPref.ItemName   = 'item'; % name of a special tag used to itemize 1D cell arrays
  DPref.StructItem = true;  % allow arrays of structs to use 'item' notation
  DPref.CellItem   = true;  % allow cell arrays to use 'item' notation
  DPref.StructTable= 'Html';
  DPref.CellTable  = 'Html';
  DPref.XmlEngine  = 'Matlab';  % use matlab provided XMLUtils
  %DPref.XmlEngine  = 'Xerces';  % use Xerces xml generator directly
  DPref.PreserveSpace = false; % Preserve or delete spaces at the beggining and the end of stings?
  RootOnly         = true;  % Input is root node only
  GlobalProcInst = [];
  GlobalComment  = [];
  GlobalDocType  = [];

  % read user preferences
  if (nargin>3)
    if (isfield(Pref, 'TableName' )),  DPref.TableName  = Pref.TableName; end
    if (isfield(Pref, 'ItemName'  )), DPref.ItemName   = Pref.ItemName;   end
    if (isfield(Pref, 'StructItem')), DPref.StructItem = Pref.StructItem; end
    if (isfield(Pref, 'CellItem'  )), DPref.CellItem   = Pref.CellItem;   end
    if (isfield(Pref, 'CellTable')),   DPref.CellTable  = Pref.CellTable; end
    if (isfield(Pref, 'StructTable')), DPref.StructTable= Pref.StructTable; end
    if (isfield(Pref, 'XmlEngine' )), DPref.XmlEngine  = Pref.XmlEngine;  end
    if (isfield(Pref, 'RootOnly'  )), RootOnly         = Pref.RootOnly;   end
    if (isfield(Pref, 'PreserveSpace')), DPref.PreserveSpace = Pref.PreserveSpace; end
  end
  if (nargin<3 || isempty(RootName)), RootName=inputname(2); end
  if (isempty(RootName)), RootName='ROOT'; end
  if (iscell(RootName)) % RootName also stores global text node data
    rName = RootName;
    RootName = char(rName{1});
    if (length(rName)>1), GlobalProcInst = char(rName{2}); end
    if (length(rName)>2), GlobalComment  = char(rName{3}); end
    if (length(rName)>3), GlobalDocType  = char(rName{4}); end
  end
  if(~RootOnly && isstruct(tree))  % if struct than deal with each field separatly
    fields = fieldnames(tree);
    for i=1:length(fields)
      field = fields{i};
      x = tree(1).(field);
      if (strcmp(field, 'COMMENT'))
        GlobalComment = x;
      elseif (strcmp(field, 'PROCESSING_INSTRUCTION'))
        GlobalProcInst = x;
      elseif (strcmp(field, 'DOCUMENT_TYPE'))
        GlobalDocType = x;
      else
        RootName = field;
        t = x;
      end
    end
    tree = t;
  end

  % Initialize jave object that will store xml data structure
  RootName = varName2str(RootName);
  if (~isempty(GlobalDocType))
    %   n = strfind(GlobalDocType, ' ');
    %   if (~isempty(n))
    %     dtype = com.mathworks.xml.XMLUtils.createDocumentType(GlobalDocType);
    %   end
    %   DOMnode = com.mathworks.xml.XMLUtils.createDocument(RootName, dtype);
    warning('xml_io_tools:write:docType', ...
      'DOCUMENT_TYPE node was encountered which is not supported yet. Ignoring.');
  end
  DOMnode = com.mathworks.xml.XMLUtils.createDocument(RootName);


  % Use recursive function to convert matlab data structure to XML
  root = DOMnode.getDocumentElement;
  struct2DOMnode(DOMnode, root, tree, DPref.ItemName, DPref);

  % Remove the only child of the root node
  root   = DOMnode.getDocumentElement;
  Child  = root.getChildNodes; % create array of children nodes
  nChild = Child.getLength;    % number of children
  if (nChild==1)
    node = root.removeChild(root.getFirstChild);
    while(node.hasChildNodes)
      root.appendChild(node.removeChild(node.getFirstChild));
    end
    while(node.hasAttributes)            % copy all attributes
      root.setAttributeNode(node.removeAttributeNode(node.getAttributes.item(0)));
    end
  end

  % Save exotic Global nodes
  if (~isempty(GlobalComment))
    DOMnode.insertBefore(DOMnode.createComment(GlobalComment), DOMnode.getFirstChild());
  end
  if (~isempty(GlobalProcInst))
    n = strfind(GlobalProcInst, ' ');
    if (~isempty(n))
      proc = DOMnode.createProcessingInstruction(GlobalProcInst(1:(n(1)-1)),...
        GlobalProcInst((n(1)+1):end));
      DOMnode.insertBefore(proc, DOMnode.getFirstChild());
    end
  end
  % Not supported yet as the code below does not work
  % if (~isempty(GlobalDocType))
  %   n = strfind(GlobalDocType, ' ');
  %   if (~isempty(n))
  %     dtype = DOMnode.createDocumentType(GlobalDocType);
  %     DOMnode.insertBefore(dtype, DOMnode.getFirstChild());
  %   end
  % end

  % save java DOM tree to XML file
  if (~isempty(filename))
    if (strcmpi(DPref.XmlEngine, 'Xerces'))
      xmlwrite_xerces(filename, DOMnode);
    else
      xmlwrite(filename, DOMnode);
    end
  end
end

  % =======================================================================
  %  === struct2DOMnode Function ===========================================
  %  =======================================================================
  function [] = struct2DOMnode(xml, parent, s, TagName, Pref)
  % struct2DOMnode is a recursive function that converts matlab's structs to
  % DOM nodes.
  % INPUTS:
  %  xml - jave object that will store xml data structure
  %  parent - parent DOM Element
  %  s - Matlab data structure to save
  %  TagName - name to be used in xml tags describing 's'
  %  Pref - preferenced
  % OUTPUT:
  %  parent - modified 'parent'

    % perform some conversions
    if (ischar(s) && min(size(s))>1) % if 2D array of characters
      s=cellstr(s);                  % than convert to cell array
    end
    % if (strcmp(TagName, 'CONTENT'))
    %   while (iscell(s) && length(s)==1), s = s{1}; end % unwrap cell arrays of length 1
    % end
    TagName  = varName2str(TagName);

    % == node is a 2D cell array ==
    % convert to some other format prior to further processing
    nDim = nnz(size(s)>1);  % is it a scalar, vector, 2D array, 3D cube, etc?
    if (iscell(s) && nDim==2 && strcmpi(Pref.CellTable, 'Matlab'))
      s = var2str(s, Pref.PreserveSpace);
    end
    if (nDim==2 && (iscell  (s) && strcmpi(Pref.CellTable,   'Vector')) || ...
                   (isstruct(s) && strcmpi(Pref.StructTable, 'Vector')))
      s = s(:);
    end
    if (nDim>2), s = s(:); end % can not handle this case well
    nItem = numel(s);
    nDim  = nnz(size(s)>1);  % is it a scalar, vector, 2D array, 3D cube, etc?

    % == node is a cell ==
    if (iscell(s)) % if this is a cell or cell array
      if ((nDim==2 && strcmpi(Pref.CellTable,'Html')) || (nDim< 2 && Pref.CellItem))
        % if 2D array of cells than can use HTML-like notation or if 1D array
        % than can use item notation
        if (strcmp(TagName, 'CONTENT')) % CONTENT nodes already have <TagName> ... </TagName>
          array2DOMnode(xml, parent, s, Pref.ItemName, Pref ); % recursive call
        else
          node = xml.createElement(TagName);   % <TagName> ... </TagName>
          array2DOMnode(xml, node, s, Pref.ItemName, Pref ); % recursive call
          parent.appendChild(node);
        end
      else % use  <TagName>...<\TagName> <TagName>...<\TagName> notation
        array2DOMnode(xml, parent, s, TagName, Pref ); % recursive call
      end
    % == node is a struct ==
    elseif (isstruct(s))  % if struct than deal with each field separatly
      if ((nDim==2 && strcmpi(Pref.StructTable,'Html')) || (nItem>1 && Pref.StructItem))
        % if 2D array of structs than can use HTML-like notation or
        % if 1D array of structs than can use 'items' notation
        node = xml.createElement(TagName);
        array2DOMnode(xml, node, s, Pref.ItemName, Pref ); % recursive call
        parent.appendChild(node);
      elseif (nItem>1) % use  <TagName>...<\TagName> <TagName>...<\TagName> notation
        array2DOMnode(xml, parent, s, TagName, Pref ); % recursive call
      else % otherwise save each struct separatelly
        fields = fieldnames(s);
        node = xml.createElement(TagName);
        for i=1:length(fields) % add field by field to the node
          field = fields{i};
          x = s.(field);
          switch field
            case {'COMMENT', 'CDATA_SECTION', 'PROCESSING_INSTRUCTION'}
              if iscellstr(x)  % cell array of strings -> add them one by one
                array2DOMnode(xml, node, x(:), field, Pref ); % recursive call will modify 'node'
              elseif ischar(x) % single string -> add it
                struct2DOMnode(xml, node, x, field, Pref ); % recursive call will modify 'node'
              else % not a string - Ignore
                warning('xml_io_tools:write:badSpecialNode', ...
                 ['Struct field named ',field,' encountered which was not a string. Ignoring.']);
              end
            case 'ATTRIBUTE' % set attributes of the node
              if (isempty(x)), continue; end
              if (isstruct(x))
                attName = fieldnames(x);       % get names of all the attributes
                for k=1:length(attName)        % attach them to the node
                  att = xml.createAttribute(varName2str(attName(k)));
                  att.setValue(var2str(x.(attName{k}),Pref.PreserveSpace));
                  node.setAttributeNode(att);
                end
              else
                warning('xml_io_tools:write:badAttribute', ...
                  'Struct field named ATTRIBUTE encountered which was not a struct. Ignoring.');
              end
            otherwise                            % set children of the node
              if ~isempty(s.(field))
                struct2DOMnode(xml, node, x, field, Pref ); % recursive call will modify 'node'
              end
          end
        end  % end for i=1:nFields
        parent.appendChild(node);
      end
    % == node is a leaf node ==
    else  % if not a struct and not a cell than it is a leaf node
      switch TagName % different processing depending on desired type of the node
        case 'COMMENT'   % create comment node
          com = xml.createComment(s);
          parent.appendChild(com);
        case 'CDATA_SECTION' % create CDATA Section
          cdt = xml.createCDATASection(s);
          parent.appendChild(cdt);
        case 'PROCESSING_INSTRUCTION' % set attributes of the node
          OK = false;
          if (ischar(s))
            n = strfind(s, ' ');
            if (~isempty(n))
              proc = xml.createProcessingInstruction(s(1:(n(1)-1)),s((n(1)+1):end));
              parent.insertBefore(proc, parent.getFirstChild());
              OK = true;
            end
          end
          if (~OK)
            warning('xml_io_tools:write:badProcInst', ...
              ['Struct field named PROCESSING_INSTRUCTION need to be',...
              ' a string, for example: xml-stylesheet type="text/css" ', ...
              'href="myStyleSheet.css". Ignoring.']);
          end
        case 'CONTENT' % this is text part of already existing node
          txt  = xml.createTextNode(var2str(s, Pref.PreserveSpace)); % convert to text
          parent.appendChild(txt);
        otherwise      % I guess it is a regular text leaf node
          txt  = xml.createTextNode(var2str(s, Pref.PreserveSpace));
          node = xml.createElement(TagName);
          node.appendChild(txt);
          parent.appendChild(node);
      end
    end % of struct2DOMnode function
  end
  % =======================================================================
  %  === array2DOMnode Function ============================================
  %  =======================================================================
  function [] = array2DOMnode(xml, parent, s, TagName, Pref)
    % Deal with 1D and 2D arrays of cell or struct. Will modify 'parent'.
    nDim = nnz(size(s)>1);  % is it a scalar, vector, 2D array, 3D cube, etc?
    switch nDim
      case 2 % 2D array
        for r=1:size(s,1)
          subnode = xml.createElement(Pref.TableName{1});
          for c=1:size(s,2)
            v = s(r,c);
            if iscell(v), v = v{1}; end
            struct2DOMnode(xml, subnode, v, Pref.TableName{2}, Pref ); % recursive call
          end
          parent.appendChild(subnode);
        end
      case 1 %1D array
        for iItem=1:numel(s)
          v = s(iItem);
          if iscell(v), v = v{1}; end
          struct2DOMnode(xml, parent, v, TagName, Pref ); % recursive call
        end
      case 0 % scalar -> this case should never be called
        if ~isempty(s) 
          if iscell(s), s = s{1}; end
          struct2DOMnode(xml, parent, s, TagName, Pref );
        end
    end
  end
  % =======================================================================
  %  === var2str Function ==================================================
  %  =======================================================================
  function str = var2str(object, PreserveSpace)
    % convert matlab variables to a string
    switch (1)
      case isempty(object)
        str = '';
      case (isnumeric(object) || islogical(object))
        if ndims(object)>2, object=object(:); end  % can't handle arrays with dimention > 2
        str=mat2str(object);           % convert matrix to a string
        % mark logical scalars with [] (logical arrays already have them) so the xml_read
        % recognizes them as MATLAB objects instead of strings. Same with sparse
        % matrices
        if ((islogical(object) && isscalar(object)) || issparse(object)),
          str = ['[' str ']'];
        end
        if (isinteger(object)),
          str = ['[', class(object), '(', str ')]'];
        end
      case iscell(object)
        if ndims(object)>2, object=object(:); end  % can't handle cell arrays with dimention > 2
        [nr nc] = size(object);
        obj2 = object;
        for i=1:length(object(:))
          str = var2str(object{i}, PreserveSpace);
          if (ischar(object{i})), object{i} = ['''' object{i} '''']; else object{i}=str; end
          obj2{i} = [object{i} ','];
        end
        for r = 1:nr, obj2{r,nc} = [object{r,nc} ';']; end
        obj2 = obj2.';
        str = ['{' obj2{:} '}'];
      case isstruct(object)
        str='';
        warning('xml_io_tools:write:var2str', ...
          'Struct was encountered where string was expected. Ignoring.');
      case isa(object, 'function_handle')
        str = ['[@' char(object) ']'];
      case ischar(object)
        str = object;
      otherwise
        str = char(object);
    end

    % string clean-up
    str=str(:); str=str.';            % make sure this is a row vector of char's
    if (~isempty(str))
      str(str<32|str==127)=' ';       % convert no-printable characters to spaces
      if (~PreserveSpace)
        str = strtrim(str);             % remove spaces from begining and the end
        str = regexprep(str,'\s+',' '); % remove multiple spaces
      end
    end
  end
  % =======================================================================
  %  === var2Namestr Function ==============================================
  %  =======================================================================
  function str = varName2str(str)
    % convert matlab variable names to a sting
    str = char(str);
    p   = strfind(str,'0x');
    if (~isempty(p))
      for i=1:length(p)
        before = str( p(i)+(0:3) );          % string to replace
        after  = char(hex2dec(before(3:4))); % string to replace with
        str = regexprep(str,before,after, 'once', 'ignorecase');
        p=p-3; % since 4 characters were replaced with one - compensate
      end
    end
    str = regexprep(str,'_COLON_',':', 'once', 'ignorecase');
    str = regexprep(str,'_DASH_' ,'-', 'once', 'ignorecase');
  end

  
  
function varargout=xmlwrite_xerces(varargin)
%XMLWRITE_XERCES Serialize an XML Document Object Model node using Xerces parser. 
%  xmlwrite_xerces(FILENAME,DOMNODE) serializes the DOMNODE to file FILENAME.
%
% The function xmlwrite_xerces is very similar the Matlab function xmlwrite 
% but works directly with the XERCES java classes (written by Apache XML 
% Project) instead of the XMLUtils class created by Mathworks. Xerces files
% are provided in standard MATLAB instalation and live in root\java\jarext
% directory. 
%
% Written by A.Amaro (02-22-2007) and generously donated to xml_io_tools. 
% This function is needed as a work-around for a bug in XMLUtils library
% which can not write CDATA SECTION nodes correctly. Also Xerces and 
% XMLUtils libraries handle namespaces differently.  
%
% Examples:
%   % See xmlwrite examples this function have almost identical behavior.
%  
% Advanced use:
%  FILENAME can also be a URN, java.io.OutputStream or java.io.Writer object
%  SOURCE can also be a SAX InputSource, JAXP Source, InputStream, or 
%    Reader object

  returnString = false;
  if length(varargin)==1
      returnString = true;
      result = java.io.StringWriter;
      source = varargin{1};
  else
      result = varargin{1};
      if ischar(result)
        % Using the XERCES classes directly, is not needed to modify the
        % filename string. So I have commented this next line
        %  result = F_xmlstringinput(result,false);
      end

      source = varargin{2};
      if ischar(source)
          source = F_xmlstringinput(source,true);
      end
  end

  % SERIALIZATION OF THE DOM DOCUMENT USING XERCES CLASSES DIRECTLY

  % 1) create the output format according to the document definitions
  % and type
  objOutputFormat = org.apache.xml.serialize.OutputFormat(source);
  set(objOutputFormat,'Indenting','on');

  % 2) create the output stream. In this case: an XML file
  objFile = java.io.File(result);
  objOutputStream = java.io.FileOutputStream(objFile);

  % 3) Create the Xerces Serializer object
  objSerializer= org.apache.xml.serialize.XMLSerializer(objOutputStream,objOutputFormat);

  % 4) Serialize to the XML files
  javaMethod('serialize',objSerializer,source);

  % 5) IMPORTANT! Delete the objects to liberate the XML file created
  objOutputStream.close;

  if returnString
      varargout{1}=char(result.toString);
  end

  % ========================================================================
  function out = F_xmlstringinput(xString,isFullSearch,varargin)
    % The function F_xmlstringinput is a copy of the private function:
    % 'xmlstringinput' that the original xmlwrite function uses.

    if isempty(xString)
        error('Filename is empty');
    elseif ~isempty(strfind(xString,'://'))
        %xString is already a URL, most likely prefaced by file:// or http://
        out = xString;
        return;
    end

    xPath=fileparts(xString);
    if isempty(xPath)
        if nargin<2 || isFullSearch
            out = which(xString);
            if isempty(out)
                error('xml:FileNotFound','File %s not found',xString);
            end
        else
            out = fullfile(pwd,xString);
        end
    else
        out = xString;
        if (nargin<2 || isFullSearch) && ~exist(xString,'file')
            %search to see if xString exists when isFullSearch
            error('xml:FileNotFound','File %s not found',xString);
        end
    end

    %Return as a URN
    if strncmp(out,'\\',2)
        % SAXON UNC filepaths need to look like file:///\\\server-name\
        out = ['file:///\',out];
    elseif strncmp(out,'/',1)
        % SAXON UNIX filepaths need to look like file:///root/dir/dir
        out = ['file://',out];
    else
        % DOS filepaths need to look like file:///d:/foo/bar
        out = ['file:///',strrep(out,'\','/')];
    end
  end
end

  
%{
function y = base64encode(x, alg, isChunked, url_safe)
%BASE64ENCODE Perform base64 encoding on a string.
% INPUT:
%   x    - block of data to be encoded.  Can be a string or a numeric
%          vector containing integers in the range 0-255.
%   alg  - Algorithm to use: can take values 'java' or 'matlab'. Optional
%          variable defaulting to 'java' which is a little faster. If
%          'java' is chosen than core of the code is performed by a call to
%          a java library. Optionally all operations can be performed using
%          matleb code.
%   isChunked - encode output into 76 character blocks. The returned 
%          encoded string is broken into lines of no more than
%          76 characters each, and each line will end with EOL. Notice that
%          if resulting string is saved as part of an xml file, those EOL's
%          are often stripped by xmlwrite funtrion prior to saving.      
%   url_safe - use Modified Base64 for URL applications ('base64url'
%   encoding) "Base64 alphabet" ([A-Za-z0-9-_=]). 
%
%
% OUTPUT:
%   y    - character array using only "Base64 alphabet" characters 
%
%   This function may be used to encode strings into the Base64 encoding
%   specified in RFC 2045 - MIME (Multipurpose Internet Mail Extensions).  
%   The Base64 encoding is designed to represent arbitrary sequences of 
%   octets in a form that need not be humanly readable.  A 65-character 
%   subset ([A-Za-z0-9+/=]) of US-ASCII is used, enabling 6 bits to be 
%   represented per printable character.
%
%   See also BASE64DECODE.
%
%   Written by Jarek Tuszynski, SAIC, jaroslaw.w.tuszynski_at_saic.com
%
%   Matlab version based on 2004 code by Peter J. Acklam
%   E-mail:      pjacklam@online.no
%   URL:         http://home.online.no/~pjacklam
%   http://home.online.no/~pjacklam/matlab/software/util/datautil/base64encode.m

  if nargin<2, alg='java';      end
  if nargin<3, isChunked=false; end
  if ~islogical(isChunked) 
    if isnumeric(isChunked)
      isChunked=(isChunked>0);
    else
      isChunked=false;
    end 
  end
  if nargin<4, url_safe=false; end
  if ~islogical(url_safe) 
    if isnumeric(url_safe)
      url_safe=(url_safe>0);
    else
      url_safe=false;
    end 
  end


  % if x happen to be a filename than read the file
  if (numel(x)<256)
    if (exist(x, 'file')==2)
      fid = fopen(x,'rb');
      x = fread(fid, 'uint8');             % read image file as a raw binary
      fclose(fid);
    end
  end

  % Perform conversion
  switch (alg)
    case 'java'
      base64 = org.apache.commons.codec.binary.Base64;
      y = base64.encodeBase64(x, isChunked); 
      if url_safe
        y = strrep(y,'=','-');
        y = strrep(y,'/','_');
      end

    case 'matlab'

      % add padding if necessary, to make the length of x a multiple of 3
      x   = uint8(x(:));
      ndbytes = length(x);                 % number of decoded bytes
      nchunks = ceil(ndbytes / 3);         % number of chunks/groups
      if rem(ndbytes, 3)>0
        x(end+1 : 3*nchunks) = 0;          % add padding
      end
      x = reshape(x, [3, nchunks]);        % reshape the data
      y = repmat(uint8(0), 4, nchunks);    % for the encoded data

      % Split up every 3 bytes into 4 pieces
      %    aaaaaabb bbbbcccc ccdddddd
      % to form
      %    00aaaaaa 00bbbbbb 00cccccc 00dddddd
      y(1,:) = bitshift(x(1,:), -2);                  % 6 highest bits of x(1,:)
      y(2,:) = bitshift(bitand(x(1,:), 3), 4);        % 2 lowest  bits of x(1,:)
      y(2,:) = bitor(y(2,:), bitshift(x(2,:), -4));   % 4 highest bits of x(2,:)
      y(3,:) = bitshift(bitand(x(2,:), 15), 2);       % 4 lowest  bits of x(2,:)
      y(3,:) = bitor(y(3,:), bitshift(x(3,:), -6));   % 2 highest bits of x(3,:)
      y(4,:) = bitand(x(3,:), 63);                    % 6 lowest  bits of x(3,:)

      % Perform the mapping
      %   0  - 25  ->  A-Z
      %   26 - 51  ->  a-z
      %   52 - 61  ->  0-9
      %   62       ->  +
      %   63       ->  /
      map = ['A':'Z', 'a':'z', '0':'9', '+/'];
      if (url_safe), map(63:64)='-_'; end
      y = map(y(:)+1);

      % Add padding if necessary.
      npbytes = 3 * nchunks - ndbytes;    % number of padding bytes
      if npbytes>0
        y(end-npbytes+1 : end) = '=';     % '=' is used for padding
      end

      % break into lines with length LineLength
      if (isChunked)
        eol = sprintf('\n');
        nebytes = numel(y);
        nlines  = ceil(nebytes / 76);     % number of lines
        neolbytes = length(eol);          % number of bytes in eol string

        % pad data so it becomes a multiple of 76 elements
        y(nebytes + 1 : 76 * nlines) = 0;
        y = reshape(y, 76, nlines);

        % insert eol strings
        y(end + 1 : end + neolbytes, :) = eol(:, ones(1, nlines));

        % remove padding, but keep the last eol string
        m = nebytes + neolbytes * (nlines - 1);
        n = (76+neolbytes)*nlines - neolbytes;
        y(m+1 : n) = [];
      end
  end

  % reshape to a row vector and make it a character array
  y = char(reshape(y, 1, numel(y)));
end
function y = base64decode(x, outfname, alg)
  %BASE64DECODE Perform base64 decoding on a string.
  %
  % INPUT:
  %   x    - block of data to be decoded.  Can be a string or a numeric  
  %          vector containing integers in the range 0-255. Any character
  %          not part of the 65-character base64 subset set is silently
  %          ignored.  Characters occuring after a '=' padding character are 
  %          never decoded. If the length of the string to decode (after 
  %          ignoring non-base64 chars) is not a multiple of 4, then a 
  %          warning is generated.
  %
  %   outfname - if provided the binary date from decoded string will be
  %          saved into a file. Since Base64 coding is often used to embbed
  %          binary data in xml files, this option can be used to extract and
  %          save them.
  %
  %   alg  - Algorithm to use: can take values 'java' or 'matlab'. Optional
  %          variable defaulting to 'java' which is a little faster. If 
  %          'java' is chosen than core of the code is performed by a call to
  %          a java library. Optionally all operations can be performed using
  %          matleb code. 
  %
  % OUTPUT:
  %   y    - array of binary data returned as uint8 
  %
  %   This function is used to decode strings from the Base64 encoding specified
  %   in RFC 2045 - MIME (Multipurpose Internet Mail Extensions).  The Base64
  %   encoding is designed to represent arbitrary sequences of octets in a form
  %   that need not be humanly readable.  A 65-character subset ([A-Za-z0-9+/=])
  %   of US-ASCII is used, enabling 6 bits to be represented per printable
  %   character.
  %
  %   See also BASE64ENCODE.
  %
  %   Written by Jarek Tuszynski, SAIC, jaroslaw.w.tuszynski_at_saic.com
  %
  %   Matlab version based on 2004 code by Peter J. Acklam
  %   E-mail:      pjacklam@online.no
  %   URL:         http://home.online.no/~pjacklam
  %   http://home.online.no/~pjacklam/matlab/software/util/datautil/base64encode.m

  if nargin<3, alg='java';  end
  if nargin<2, outfname=''; end

  % if x happen to be a filename than read the file
  if (numel(x)<256)
    if (exist(x, 'file')==2)
      fid = fopen(x,'rb');
      x = fread(fid, 'uint8');   
      fclose(fid);
    end
  end
  x = uint8(x(:)); % unify format

  % Perform conversion
  switch (alg)
    case 'java' 
      base64 = org.apache.commons.codec.binary.Base64;
      y = base64.decode(x);
      y = mod(int16(y),256); % convert from int8 to uint8
    case 'matlab'
      %  Perform the mapping
      %   A-Z  ->  0  - 25
      %   a-z  ->  26 - 51
      %   0-9  ->  52 - 61
      %   + -  ->  62       '-' is URL_SAFE alternative
      %   / _  ->  63       '_' is URL_SAFE alternative
      map = uint8(zeros(1,256)+65);
      map(uint8(['A':'Z', 'a':'z', '0':'9', '+/=']))= 0:64;
      map(uint8('-_'))= 62:63;  % URL_SAFE alternatives
      x = map(x);  % mapping

      x(x>64)=[]; % remove non-base64 chars
      if rem(numel(x), 4)
        warning('Length of base64 data not a multiple of 4; padding input.');
      end
      x(x==64)=[]; % remove padding characters

      % add padding and reshape
      nebytes = length(x);         % number of encoded bytes
      nchunks = ceil(nebytes/4);   % number of chunks/groups
      if rem(nebytes, 4)>0
        x(end+1 : 4*nchunks) = 0;  % add padding
      end
      x = reshape(uint8(x), 4, nchunks);
      y = repmat(uint8(0), 3, nchunks);            % for the decoded data

      % Rearrange every 4 bytes into 3 bytes
      %    00aaaaaa 00bbbbbb 00cccccc 00dddddd
      % to form
      %    aaaaaabb bbbbcccc ccdddddd
      y(1,:) = bitshift(x(1,:), 2);                 % 6 highest bits of y(1,:)
      y(1,:) = bitor(y(1,:), bitshift(x(2,:), -4)); % 2 lowest  bits of y(1,:)
      y(2,:) = bitshift(x(2,:), 4);                 % 4 highest bits of y(2,:)
      y(2,:) = bitor(y(2,:), bitshift(x(3,:), -2)); % 4 lowest  bits of y(2,:)
      y(3,:) = bitshift(x(3,:), 6);                 % 2 highest bits of y(3,:)
      y(3,:) = bitor(y(3,:), x(4,:));               % 6 lowest  bits of y(3,:)

      % remove extra padding
      switch rem(nebytes, 4)
        case 2
          y = y(1:end-2);
        case 3
          y = y(1:end-1);
      end
  end

  % reshape to a row vector and make it a character array
  y = uint8(reshape(y, 1, numel(y)));

  % save to file if needed
  if ~isempty(outfname)
    fid = fopen(outfname,'wb');
    fwrite(fid, y, 'uint8');  
    fclose(fid);
  end
end
%}



%{
This works only for special structures... grummel. 

% struct2xml - function:
% __________________________________________________________________________________________________

function varargout = struct2xml( s, varargin )
%Convert a MATLAB structure into a xml file 
% [ ] = struct2xml( s, file )
% xml = struct2xml( s )
%
% A structure containing:
% s.XMLname.Attributes.attrib1 = "Some value";
% s.XMLname.Element.Text = "Some text";
% s.XMLname.DifferentElement{1}.Attributes.attrib2 = "2";
% s.XMLname.DifferentElement{1}.Text = "Some more text";
% s.XMLname.DifferentElement{2}.Attributes.attrib3 = "2";
% s.XMLname.DifferentElement{2}.Attributes.attrib4 = "1";
% s.XMLname.DifferentElement{2}.Text = "Even more text";
%
% Will produce:
% <XMLname attrib1="Some value">
%   <Element>Some text</Element>
%   <DifferentElement attrib2="2">Some more text</Element>
%   <DifferentElement attrib3="2" attrib4="1">Even more text</DifferentElement>
% </XMLname>
%
% Please note that the following strings are substituted
% '_dash_' by '-', '_colon_' by ':' and '_dot_' by '.'
%
% Written by W. Falkena, ASTI, TUDelft, 27-08-2010
% On-screen output functionality added by P. Orth, 01-12-2010
% Multiple space to single space conversion adapted for speed by T. Lohuis, 11-04-2011
% Val2str subfunction bugfix by H. Gsenger, 19-9-2011
    
    if (nargin ~= 2)
        if(nargout ~= 1 || nargin ~= 1)
            error(['Supported function calls:' sprintf('\n')...
                   '[ ] = struct2xml( s, file )' sprintf('\n')...
                   'xml = struct2xml( s )']);
        end
    end

    if(nargin == 2)
        file = varargin{1};

        if (isempty(file))
            error('Filename can not be empty');
        end

        if (isempty(strfind(file,'.xml')))
            file = [file '.xml'];
        end
    end
    
    if (~isstruct(s))
        error([inputname(1) ' is not a structure']);
    end
    
    if (length(fieldnames(s)) > 1)
        error(['Error processing the structure:' sprintf('\n') 'There should be a single field in the main structure.']);
    end
    xmlname = fieldnames(s);
    xmlname = xmlname{1};
    
    %substitute special characters
    xmlname_sc = xmlname;
    xmlname_sc = strrep(xmlname_sc,'_dash_','-');
    xmlname_sc = strrep(xmlname_sc,'_colon_',':');
    xmlname_sc = strrep(xmlname_sc,'_dot_','.');

    %create xml structure
    docNode = com.mathworks.xml.XMLUtils.createDocument(xmlname_sc);

    %process the rootnode
    docRootNode = docNode.getDocumentElement;

    %append childs
    parseStruct(s.(xmlname),docNode,docRootNode,[inputname(1) '.' xmlname '.']);

    if(nargout == 0)
        %save xml file
        xmlwrite(file,docNode);
    else
        varargout{1} = xmlwrite(docNode);
    end  
end

% ----- Subfunction parseStruct -----
function [] = parseStruct(s,docNode,curNode,pName)
    
    fnames = fieldnames(s);
    for i = 1:length(fnames)
        curfield = fnames{i};
        
        %substitute special characters
        curfield_sc = curfield;
        curfield_sc = strrep(curfield_sc,'_dash_','-');
        curfield_sc = strrep(curfield_sc,'_colon_',':');
        curfield_sc = strrep(curfield_sc,'_dot_','.');
        
        if (strcmp(curfield,'Attributes'))
            %Attribute data
            if (isstruct(s.(curfield)))
                attr_names = fieldnames(s.Attributes);
                for a = 1:length(attr_names)
                    cur_attr = attr_names{a};
                    
                    %substitute special characters
                    cur_attr_sc = cur_attr;
                    cur_attr_sc = strrep(cur_attr_sc,'_dash_','-');
                    cur_attr_sc = strrep(cur_attr_sc,'_colon_',':');
                    cur_attr_sc = strrep(cur_attr_sc,'_dot_','.');
                    
                    [cur_str,succes] = val2str(s.Attributes.(cur_attr));
                    if (succes)
                        curNode.setAttribute(cur_attr_sc,cur_str);
                    else
                        disp(['Warning. The text in ' pName curfield '.' cur_attr ' could not be processed.']);
                    end
                end
            else
                disp(['Warning. The attributes in ' pName curfield ' could not be processed.']);
                disp(['The correct syntax is: ' pName curfield '.attribute_name = ''Some text''.']);
            end
        elseif (strcmp(curfield,'Text'))
            %Text data
            [txt,succes] = val2str(s.Text);
            if (succes)
                curNode.appendChild(docNode.createTextNode(txt));
            else
                disp(['Warning. The text in ' pName curfield ' could not be processed.']);
            end
        else
            %Sub-element
            if (isstruct(s.(curfield)))
                %single element
                curElement = docNode.createElement(curfield_sc);
                curNode.appendChild(curElement);
                parseStruct(s.(curfield),docNode,curElement,[pName curfield '.'])
            elseif (iscell(s.(curfield)))
                %multiple elements
                for c = 1:length(s.(curfield))
                    curElement = docNode.createElement(curfield_sc);
                    curNode.appendChild(curElement);
                    if (isstruct(s.(curfield){c}))
                        parseStruct(s.(curfield){c},docNode,curElement,[pName curfield '{' num2str(c) '}.'])
                    else
                        disp(['Warning. The cell ' pName curfield '{' num2str(c) '} could not be processed, since it contains no structure.']);
                    end
                end
            else
                %eventhough the fieldname is not text, the field could
                %contain text. Create a new element and use this text
                curElement = docNode.createElement(curfield_sc);
                curNode.appendChild(curElement);
                [txt,succes] = val2str(s.(curfield));
                if (succes)
                    curElement.appendChild(docNode.createTextNode(txt));
                else
                    disp(['Warning. The text in ' pName curfield ' could not be processed.']);
                end
            end
        end
    end
end

%----- Subfunction val2str -----
function [str,succes] = val2str(val)
    
    succes = true;
    str = [];
    
    if (isempty(val))
        return; %bugfix from H. Gsenger
    elseif (ischar(val))
        %do nothing
    elseif (isnumeric(val))
        val = num2str(val);
    else
        succes = false;
    end
    
    if (ischar(val))
        %add line breaks to all lines except the last (for multiline strings)
        lines = size(val,1);
        val = [val char(sprintf('\n')*[ones(lines-1,1);0])];
        
        %transpose is required since indexing (i.e., val(nonspace) or val(:)) produces a 1-D vector. 
        %This should be row based (line based) and not column based.
        valt = val';
        
        remove_multiple_white_spaces = true;
        if (remove_multiple_white_spaces)
            %remove multiple white spaces using isspace, suggestion of T. Lohuis
            whitespace = isspace(val);
            nonspace = (whitespace + [zeros(lines,1) whitespace(:,1:end-1)])~=2;
            nonspace(:,end) = [ones(lines-1,1);0]; %make sure line breaks stay intact
            str = valt(nonspace');
        else
            str = valt(:);
        end
    end
end


% xml2struct - function:
% __________________________________________________________________________________________________

function [ s ] = xml2struct( file )
%Convert xml file into a MATLAB structure
% [ s ] = xml2struct( file )
%
% A file containing:
% <XMLname attrib1="Some value">
%   <Element>Some text</Element>
%   <DifferentElement attrib2="2">Some more text</Element>
%   <DifferentElement attrib3="2" attrib4="1">Even more text</DifferentElement>
% </XMLname>
%
% Will produce:
% s.XMLname.Attributes.attrib1 = "Some value";
% s.XMLname.Element.Text = "Some text";
% s.XMLname.DifferentElement{1}.Attributes.attrib2 = "2";
% s.XMLname.DifferentElement{1}.Text = "Some more text";
% s.XMLname.DifferentElement{2}.Attributes.attrib3 = "2";
% s.XMLname.DifferentElement{2}.Attributes.attrib4 = "1";
% s.XMLname.DifferentElement{2}.Text = "Even more text";
%
% Please note that the following characters are substituted
% '-' by '_dash_', ':' by '_colon_' and '.' by '_dot_'
%
% Written by W. Falkena, ASTI, TUDelft, 21-08-2010
% Attribute parsing speed increased by 40% by A. Wanner, 14-6-2011
% Added CDATA support by I. Smirnov, 20-3-2012
%
% Modified by X. Mo, University of Wisconsin, 12-5-2012

    if (nargin < 1)
        clc;
        help xml2struct
        return
    end
    
    if isa(file, 'org.apache.xerces.dom.DeferredDocumentImpl') || isa(file, 'org.apache.xerces.dom.DeferredElementImpl')
        % input is a java xml object
        xDoc = file;
    else
        %check for existance
        if (exist(file,'file') == 0)
            %Perhaps the xml extension was omitted from the file name. Add the
            %extension and try again.
            if (isempty(strfind(file,'.xml')))
                file = [file '.xml'];
            end
            
            if (exist(file,'file') == 0)
                error(['The file ' file ' could not be found']);
            end
        end
        %read the xml file
        xDoc = xmlread(file);
    end
    
    %parse xDoc into a MATLAB structure
    s = parseChildNodes(xDoc);
    
end

% ----- Subfunction parseChildNodes -----
function [children,ptext,textflag] = parseChildNodes(theNode)
    % Recurse over node children.
    children = struct;
    ptext = struct; textflag = 'Text';
    if hasChildNodes(theNode)
        childNodes = getChildNodes(theNode);
        numChildNodes = getLength(childNodes);

        for count = 1:numChildNodes
            theChild = item(childNodes,count-1);
            [text,name,attr,childs,textflag] = getNodeData(theChild);
            
            if (~strcmp(name,'#text') && ~strcmp(name,'#comment') && ~strcmp(name,'#cdata_dash_section'))
                %XML allows the same elements to be defined multiple times,
                %put each in a different cell
                if (isfield(children,name))
                    if (~iscell(children.(name)))
                        %put existsing element into cell format
                        children.(name) = {children.(name)};
                    end
                    index = length(children.(name))+1;
                    %add new element
                    children.(name){index} = childs;
                    if(~isempty(fieldnames(text)))
                        children.(name){index} = text; 
                    end
                    if(~isempty(attr)) 
                        children.(name){index}.('Attributes') = attr; 
                    end
                else
                    %add previously unknown (new) element to the structure
                    children.(name) = childs;
                    if(~isempty(text) && ~isempty(fieldnames(text)))
                        children.(name) = text; 
                    end
                    if(~isempty(attr)) 
                        children.(name).('Attributes') = attr; 
                    end
                end
            else
                ptextflag = 'Text';
                if (strcmp(name, '#cdata_dash_section'))
                    ptextflag = 'CDATA';
                elseif (strcmp(name, '#comment'))
                    ptextflag = 'Comment';
                end
                
                %this is the text in an element (i.e., the parentNode) 
                if (~isempty(regexprep(text.(textflag),'[\s]*','')))
                    if (~isfield(ptext,ptextflag) || isempty(ptext.(ptextflag)))
                        ptext.(ptextflag) = text.(textflag);
                    else
                        %what to do when element data is as follows:
                        %<element>Text <!--Comment--> More text</element>
                        
                        %put the text in different cells:
                        % if (~iscell(ptext)) ptext = {ptext}; end
                        % ptext{length(ptext)+1} = text;
                        
                        %just append the text
                        ptext.(ptextflag) = [ptext.(ptextflag) text.(textflag)];
                    end
                end
            end
            
        end
    end
end

% ----- Subfunction getNodeData -----
function [text,name,attr,childs,textflag] = getNodeData(theNode)
    % Create structure of node info.
    
    %make sure name is allowed as structure name
    name = toCharArray(getNodeName(theNode))';
    name = strrep(name, '-', '_dash_');
    name = strrep(name, ':', '_colon_');
    name = strrep(name, '.', '_dot_');

    attr = parseAttributes(theNode);
    if (isempty(fieldnames(attr))) 
        attr = []; 
    end
    
    %parse child nodes
    [childs,text,textflag] = parseChildNodes(theNode);
    
    if (isempty(fieldnames(childs)) && isempty(fieldnames(text)))
        %get the data of any childless nodes
        % faster than if any(strcmp(methods(theNode), 'getData'))
        % no need to try-catch (?)
        % faster than text = char(getData(theNode));
        text.(textflag) = toCharArray(getTextContent(theNode))';
    end
    
end

% ----- Subfunction parseAttributes -----
function attributes = parseAttributes(theNode)
    % Create attributes structure.

    attributes = struct;
    if hasAttributes(theNode)
       theAttributes = getAttributes(theNode);
       numAttributes = getLength(theAttributes);

       for count = 1:numAttributes
            %attrib = item(theAttributes,count-1);
            %attr_name = regexprep(char(getName(attrib)),'[-:.]','_');
            %attributes.(attr_name) = char(getValue(attrib));

            %Suggestion of Adrian Wanner
            str = toCharArray(toString(item(theAttributes,count-1)))';
            k = strfind(str,'='); 
            attr_name = str(1:(k(1)-1));
            attr_name = strrep(attr_name, '-', '_dash_');
            attr_name = strrep(attr_name, ':', '_colon_');
            attr_name = strrep(attr_name, '.', '_dot_');
            attributes.(attr_name) = str((k(1)+2):(end-1));
       end
    end
end

%}


