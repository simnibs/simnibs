function [CATrel, CATver]  = cat_version(varargin)
% check for CAT revision
%
% FORMAT [CATrel, CATver] = cat_version
% FORMAT CATver = cat_version % short version
% FORMAT cat_version('[ss]fnbanner'[,str,ver]) % display banner
% 
% This function will retrieve the CAT release and version and is a
% modified version of spm('version')
%_______________________________________________________________________
% Christian Gaser
% $Id: cat_version.m 1250 2017-12-20 16:17:28Z gaser $

persistent CAT_VER;
v = CAT_VER;

if isempty(CAT_VER)
    v = struct('Name','','Version','','Release','','Date','');
    % try Contents.txt and then Contents.m
    try
        %% v1 based on the Content file
        vfile = fullfile(spm('Dir'),'toolbox','cat12','Contents.txt');
        if ~exist(vfile,'file')
          vfile = fullfile(spm('Dir'),'toolbox','cat12','Contents.m');
        end
        fid = fopen(vfile,'rt');
        if fid == -1, error('Can''t open %s.',vfile); end
        l1 = fgetl(fid); l2 = fgetl(fid);
        fclose(fid);
        l1 = strtrim(l1(2:end)); l2 = strtrim(l2(2:end));
        t  = textscan(l2,'%s','delimiter',' '); t = t{1};
        v.Name    = l1;
        v.Date    = t{4};
        v.Version = t{2};
        v.Release = t{3}(2:end-1);
        %v.User    = 'gaser';
        
        %% v2 based on the CHANGES.txt
        vfile = fullfile(spm('Dir'),'toolbox','cat12','CHANGES.txt');
        if exist(vfile,'file')
          fid = fopen(vfile,'rt'); ct = textscan(fid,'%s','delimiter','\n'); fclose(fid);
          t  = textscan(ct{1}{2},'%s','delimiter',' '); t = t{1};
          v2.Date     = t{5};
          v2.Version  = t{1}(2:end); 
          %v2.User     = t{3};

          if str2double(v2.Version)>str2double(v.Version), v=cat_io_updateStruct(v,v2); end
        end
    catch
        error('Can''t obtain CAT Revision information.');
    end
    CAT_VER = v;
end

if nargin>0, Action = varargin{1}; else Action = ''; end
switch Action
  case {'fnbanner','sfnbanner','ssfnbanner'}  %-Text banners for functions
    %=======================================================================
    % SPMid = spm('FnBanner', Fn,FnV)
    % SPMid = spm('SFnBanner',Fn,FnV)
    % SPMid = spm('SSFnBanner',Fn,FnV)
    %-----------------------------------------------------------------------
    time = spm('time');
    str  = sprintf('%s(%s)',v.Release,v.Version); %spm('Ver');
    if nargin>=2, str = [str,': ',varargin{2}]; end
    if nargin>=3 
        v = regexp(varargin{3},'\$Rev: (\d*) \$','tokens','once');
        if ~isempty(v)
            str = [str,' (v',v{1},')'];
        else
            str = [str,' (v',varargin{3},')'];
        end
    end

    switch lower(Action)
    case 'fnbanner'
        tab = '';
        wid = 72;
        lch = '=';
    case 'sfnbanner'
        tab = sprintf('\t');
        wid = 72-8;
        lch = '-';
    case 'ssfnbanner'
        tab = sprintf('\t\t');
        wid = 72-2*8;
        lch = '-';
    end

    fprintf('\n%s%s',tab,str)
    fprintf('%s',repmat(' ',1,wid-length([str,time])))
    fprintf('%s\n%s',time,tab)
    fprintf('%s',repmat(lch,1,wid)),fprintf('\n')
end

CATrel = v.Release;
CATver = v.Version;
