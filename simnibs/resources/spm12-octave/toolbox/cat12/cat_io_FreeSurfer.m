function varargout=cat_io_FreeSurfer(action,varargin)
% ______________________________________________________________________
% 
% Read/Write FreeSurfer Data. 
%
% Use FreeSurfer in/output functions created by Bruce Fischl, Doug Greve,
% Thomas Yeo, and Mert Sabuncu.
%
%   varargout = cat_io_FreeSurfer(action,varargin)
%
% * surface meshs:
%   cat_io_FreeSurfer('write_surf',fname,vertices,faces);
%   S = cat_io_FreeSurfer('read_surf',fname);
%
% * surface data:
%   cat_io_FreeSurfer('write_surf_data',fname,cdata);
%   cdata = cat_io_FreeSurfer('read_surf_data',fname);
%
% * surface atlases:
%   [vertices, label, colortable] = 
%      cat_io_FreeSurfer('read_annotation',fname);
%   cat_io_FreeSurfer('write_annotation', ...
%      fname, vertices, label, colortable);
%
% * GIFTI to FreeSurfer / FreeSurfer to GIFTI: 
%   [P] = cat_io_FreeSurfer('gii2fs',fname);
%   [P] = cat_io_FreeSurfer('gii2fs',...
%           struct('data',{fnames},'delete',[0|1]));
%
%   [P] = cat_io_FreeSurfer('fs2gii',fname);
%   [P] = cat_io_FreeSurfer('fs2gii',...
%           struct('data',{fnames},'cdata',{cfnames},'delete',[0|1]));
% ______________________________________________________________________
%
%   Robert Dahnke (robert.dahnke@uni-jena.de)
%   Structural Brain Mapping Group (http://dbm.neuro.uni-jena.de/)
%   Department of Neurology
%   University Jena
% ______________________________________________________________________
% $Id: cat_io_FreeSurfer.m 1233 2017-12-03 23:26:25Z gaser $ 

  if ~exist('action','var'), help cat_io_FreeSurfer; return; end
  if nargin==0, varargin{1} = struct(); end
  
  switch action
    case {'fs2gii','read_surf','read_surf_data'}
      [pp,ff,ee] = spm_fileparts(varargin{1}); 
      if any(strcmp({'.gii','.nii','.mat','.m'},ee))
        error('cat_io_FreeSurfer:WrongInput','FreeSurfer file format only! No filenames GIFTI (*.gii), NIFTI (*.nii), or Matlab (*.mat,*.m) are allowed! ');
      end
  end
  
  switch action
    %case 'FSatlas2cat'
    %  varargout{1} = cat_surf_FSannotation2CAT(varargin{1});
    case 'gii2fs'
      if nargout==1
        varargout{1} = gii2fs(varargin{1}); 
      elseif nargout==2
        [varargout{1},varargout{2}] = gii2fs(varargin{1}); 
      else
        gii2fs(varargin{1})  
      end
    case 'fs2gii'
      if nargout>0
        varargout{1} = fs2gii(varargin{1}); 
      else
        fs2gii(varargin{1})
      end
    case 'read_annotation'
      if nargin==2
        [varargout{1}, varargout{2}, varargout{3}] = Read_Brain_Annotation(varargin{1}); 
      else
        [varargout{1}, varargout{2}, varargout{3}] = Read_Brain_Annotation(varargin{1}, varargin(2:end)); 
      end
      varargout{4} = [{'ROIid'},{'ROIname'};num2cell(varargout{3}.table(1:end,5)),varargout{3}.struct_names];
    case 'write_annotation' 
      
      %write_annotation(filename, vertices, label, ct)
      write_annotation(varargin{1}, varargin{2}, varargin{3}, varargin{4})
    case 'write_surf'
      write_surf(varargin{1}, varargin{2}.vertices, varargin{2}.faces);
    case 'read_surf'
      [varargout{1}.vertices,varargout{1}.faces] = read_surf(varargin{1}); 
      varargout{1}.faces = varargout{1}.faces+1;   
    case 'write_surf_data'
      write_curv(varargin{1},varargin{2});
    case 'read_surf_data'
      [varargout{1},varargout{2}] = read_curv(varargin{1});
    otherwise
      error(['cat_io_FreeSurfer:unknownAction','Unknown action ''%s''!\n' ...
             'Use ''write_surf'',''read_surf'',''write_surf_data'',''read_surf_data'',''gii2fs'',''fs2gii'',''read_annotation'',''write_annotation''.\n'],action);
  end

end

function job = getjob(job0,sel)
  if isstruct(job0)
    job = job0; 
  else 
    job.data = job0;
  end
  if ~isfield(job,'data') || isempty(job.data)
    job.data = spm_select(inf,'any','Select surface','','',sel);
  end
  if isempty(job.data), return; end
  
  job.data  = cellstr(job.data);
  if isfield(job,'cdata'), job.cdata = cellstr(job.cdata); end
  if isfield(job,'cdata') && isfield(job,'data') && ...
      numel(job.cdata) ~= numel(job.data)
    error('cat_io_FreeSurfer:getjob:data','Number of surface meshes and textures have to be equivalent'); 
  end
  
  for si=1:numel(job.data)
    if isfield(job,'cdata')
      [pp,ff,ee] = spm_fileparts(job.cdata{si});
    else
      [pp,ff,ee] = spm_fileparts(job.data{si});
    end
    def.fname{si} = strrep(fullfile(pp,[ff ee]),'.gii',''); 
  end
  
  def.verb    = 0; 
  def.delete  = 0; 
  def.merge   = 0;
  
  job = cat_io_checkinopt(job,def);
end
function varargout = gii2fs(varargin)
% convert gifti surfaces to FreeSurfer 
  job = getjob(varargin,'[lr]h.*.gii');
  
  surfname = cell(numel(job.data),1); 
  curfname = cell(numel(job.data),1); 
  for si=1:numel(job.data)
    [pp,ff] = spm_fileparts(job.data{si});
    sinfo  = cat_surf_info(job.data{si}); 
    
    CS = gifti(job.data{si});
    if isfield(CS,'vertices') && isfield(CS,'faces')
      switch sinfo.texture
        case {'sphere','central','hull','inner','outer'}
          surfname{si} = char(cat_surf_rename(sinfo,'ee',''));
        otherwise
          surfname{si} = char(cat_surf_rename(sinfo,'dataname',[sinfo.texture '_surface'],'ee',''));
      end
      write_surf(char(surfname{si}), CS.vertices , CS.faces);
    else
      surfname{si} = ''; 
    end
    if isfield(CS,'cdata')
      curfname{si} = fullfile(pp,ff);
      
      write_curv(curfname{si}, double(CS.cdata));
    else
      curfname{si} = ''; 
    end
    if job.delete
      delete(job.data{si}); 
    end
  end
  
  if nargout>0
    varargout{1} = surfname;
  end
  if nargout>1
    varargout{2} = curfname;
  end
    
end

function varargout = fs2gii(varargin)
% convert FreeSurfer surfaces meshes/data files to a gifti
  job = getjob(varargin{1},'[lr]h.*');
  
  for si=1:numel(job.data)
    [pp,ff,ee] = spm_fileparts(job.data{si}); %#ok<ASGLU>
    switch ee
      case '.gii'
        S = gifti(job.data{si});
      otherwise
        %if isfield(job,'cdata') 
          try
            [vertices,faces] = read_surf(job.data{si}); 
            S.vertices = vertices; 
            S.faces    = faces; 
          end
        %else
        %  try %#ok<TRYNC>
        %    S.cdata = read_curv(job.cdata{si});   
        %  end
        %end
    end  
    
    if isfield(job,'cdata') 
      S.cdata = read_curv(job.cdata{si});     
    end
    
    job.fname{si} = [job.fname{si} '.gii'];
    save(gifti(S),job.fname{si}); 
    
    if job.delete
      delete(job.data{si});
    end
  end
  
  if nargout>0
    varargout{1} = job.fname;
  end
end

function annots = cat_surf_FSannotation2CAT(job)
% -- in development --
% Read FreeSurfer average atlas maps and save them as texture with csv 
% and xml data.
% > convertation to SPM ROI format ...

  def.trerr     = 0; 
  def.verb      = cat_get_defaults('extopts.verb'); 
  def.debug     = cat_get_defaults('extopts.verb')>2;
  def.fsavgDir  = fullfile(spm('dir'),'toolbox','cat12','templates_surfaces');
  if ismac
    def.FSDir   = cat_vol_findfiles('/Applications','Freesurfer*',struct('depth',1,'dirs',1));
    def.FSDir   = def.FSDir{end};
  end
  def.FSsub     = cat_vol_findfiles(fullfile(def.FSDir,'subjects'),'fsaverage',struct('depth',1,'dirs',1));
  job = cat_io_checkinopt(job,def);
  
  for FSsubi=1:numel(job.FSsub)
    annots{FSsubi} = cat_vol_findfiles(fullfile(job.FSsub{FSsubi},'label'),'*.annot');
     
    for annotsi=1:numel(annots{FSsubi})
      %%
      [pp,ff,ee] = fileparts(annots{FSsubi}{annotsi}); side=ff(1:2);
      [vertices, label, colortable] = cat_io_FreeSurfer('read_annotation',annots{FSsubi}{annotsi}); 
    
      Pfsavg = fullfile(def.fsavgDir,sprintf('%s.central.freesurfer.gii',side));      % fsaverage central
      
      CS = gifti(Pfsavg);
    
    end
    
  end
  
end

function write_surf(fname, vert, face)
% write_surf - FreeSurfer I/O function to write a surface file
% 
% write_surf(fname, vert, face)
% 
% writes a surface triangulation into a binary file
% fname - name of file to write
% vert  - Nx3 matrix of vertex coordinates
% face  - Mx3 matrix of face triangulation indices
% 
% The face matrix here must be matlab compatible
% (no zero indices).  It is converted to FreeSurfer
% indices that start at zero.
% 
% See also freesurfer_read_surf, freesurfer_write_curv, freesurfer_write_wfile

  if(nargin ~= 3)
    fprintf('USAGE: freesurfer_write_surf(fname, vert, face)\n');
    return;
  end

  if size(vert,2) ~= 3,
      error('vert must be Nx3 matrix');
  end

  if size(face,2) ~= 3,
      error('face must be Mx3 matrix');
  end

  %fprintf('...subtracting 1 from face indices for FreeSurfer compatibility.\n');
  face = face - 1;

  % open it as a big-endian file
  fid = fopen(fname, 'wb', 'b') ;

  TRIANGLE_FILE_MAGIC_NUMBER = 16777214 ;
  fwrite3(fid, TRIANGLE_FILE_MAGIC_NUMBER);

  vnum = size(vert,1) ;  % number of vertices
  fnum = size(face,1) ;  % number of faces

  % Ouput a couple of text lines with creation date
  fprintf(fid,'created by %s on %s\n\n',getenv('USER'),datestr(now)); % creation date 

  fwrite(fid, vnum,'int32');
  fwrite(fid, fnum,'int32');

  % reshape vert into column array and write
  vert = reshape(vert',size(vert,1)*size(vert,2),1);
  fwrite(fid, vert,'float32');

  % reshape face into column array and write
  face = reshape(face',size(face,1)*size(face,2),1);
  fwrite(fid, face,'int32');

  fclose(fid) ;

end

function [vertex_coords, faces] = read_surf(fname)
  %
  % [vertex_coords, faces] = read_surf(fname)
  % reads a the vertex coordinates and face lists from a surface file
  % note that reading the faces from a quad file can take a very long
  % time due to the goofy format that they are stored in. If the faces
  % output variable is not specified, they will not be read so it 
  % should execute pretty quickly.
  %


  %
  % read_surf.m
  %
  % Original Author: Bruce Fischl
  % CVS Revision Info:
  %    $Author: gaser $
  %    $Date: 2017-12-04 00:26:25 +0100 (Mon, 04 Dec 2017) $
  %    $Revision: 1233 $
  %
  % Copyright ?? 2011 The General Hospital Corporation (Boston, MA) "MGH"
  %
  % Terms and conditions for use, reproduction, distribution and contribution
  % are found in the 'FreeSurfer Software License Agreement' contained
  % in the file 'LICENSE' found in the FreeSurfer distribution, and here:
  %
  % https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
  %
  % Reporting: freesurfer@nmr.mgh.harvard.edu
  %


  %fid = fopen(fname, 'r') ;
  %nvertices = fscanf(fid, '%d', 1);
  %all = fscanf(fid, '%d %f %f %f %f\n', [5, nvertices]) ;
  %curv = all(5, :)' ;

  % open it as a big-endian file


  %QUAD_FILE_MAGIC_NUMBER =  (-1 & 0x00ffffff) ;
  %NEW_QUAD_FILE_MAGIC_NUMBER =  (-3 & 0x00ffffff) ;

  TRIANGLE_FILE_MAGIC_NUMBER  =  16777214 ;
  QUAD_FILE_MAGIC_NUMBER      =  16777215 ;

  if ~exist(fname,'file')
    error('MATLAB:cat_io_FreeSurfer:read_surf','mesh file %s does not exist.', fname) ;
  end
  fid = fopen(fname, 'rb', 'b') ;
  if (fid < 0)
    error('MATLAB:cat_io_FreeSurfer:read_surf','could not open mesh file %s.', fname) ;
  end
  magic = fread3(fid) ;

  if(magic == QUAD_FILE_MAGIC_NUMBER)
    vnum = fread3(fid) ;
    fnum = fread3(fid) ;
    vertex_coords = fread(fid, vnum*3, 'int16') ./ 100 ; 
    faces = zeros(fnum,4,'single');
    if (nargout > 1)
      for i=1:fnum
        for n=1:4
          faces(i,n) = fread3(fid) ;
        end
      end
    end
  elseif (magic == TRIANGLE_FILE_MAGIC_NUMBER)
    fgets(fid) ;
    fgets(fid) ;
    vnum = fread(fid, 1, 'int32') ;
    fnum = fread(fid, 1, 'int32') ;
    vertex_coords = fread(fid, vnum*3, 'float32') ; 
    faces = fread(fid, fnum*3, 'int32') ;
    faces = reshape(faces, 3, fnum)' ;
  end
  if min(faces(:))==0, faces=faces+1; end
  vertex_coords = reshape(vertex_coords, 3, vnum)' ;
  fclose(fid) ;
end

function fwrite3(fid, val)
%
% fwrite3.m
%
% Original Author: Bruce Fischl
% CVS Revision Info:
%    $Author: gaser $
%    $Date: 2017-12-04 00:26:25 +0100 (Mon, 04 Dec 2017) $
%    $Revision: 1233 $
%
% Copyright ?? 2011 The General Hospital Corporation (Boston, MA) "MGH"
%
% Terms and conditions for use, reproduction, distribution and contribution
% are found in the 'FreeSurfer Software License Agreement' contained
% in the file 'LICENSE' found in the FreeSurfer distribution, and here:
%
% https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
%
% Reporting: freesurfer@nmr.mgh.harvard.edu
%

% write a 3 byte integer out of a file
%fwrite(fid, val, '3*uchar') ;
  
  b1 = bitand(bitshift(val, -16), 255) ;
  b2 = bitand(bitshift(val, -8), 255) ;
  b3 = bitand(val, 255) ; 
  fwrite(fid, b1, 'uchar') ;
  fwrite(fid, b2, 'uchar') ;
  fwrite(fid, b3, 'uchar') ;
end

function [retval] = fread3(fid)
  % [retval] = fd3(fid)
  % read a 3 byte integer out of a file


  %
  % fread3.m
  %
  % Original Author: Bruce Fischl
  % CVS Revision Info:
  %    $Author: gaser $
  %    $Date: 2017-12-04 00:26:25 +0100 (Mon, 04 Dec 2017) $
  %    $Revision: 1233 $
  %
  % Copyright ?? 2011 The General Hospital Corporation (Boston, MA) "MGH"
  %
  % Terms and conditions for use, reproduction, distribution and contribution
  % are found in the 'FreeSurfer Software License Agreement' contained
  % in the file 'LICENSE' found in the FreeSurfer distribution, and here:
  %
  % https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
  %
  % Reporting: freesurfer@nmr.mgh.harvard.edu
  %

  b1 = fread(fid, 1, 'uchar') ;
  b2 = fread(fid, 1, 'uchar') ;
  b3 = fread(fid, 1, 'uchar') ;
  retval = bitshift(b1, 16) + bitshift(b2,8) + b3 ;

end

function err = write_wfile(fname, w, v)
% err = write_wfile(fname, w, <v>)
% 
% writes a vector into a binary 'w' file
%  fname - name of file to write to
%  w     - vector of values to be written
%  v     - 0-based vertex numbers 
%          (assumes 0 to N-1 if not present or empty).
%
% See also read_wfile.
%


%
% write_wfile.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: gaser $
%    $Date: 2017-12-04 00:26:25 +0100 (Mon, 04 Dec 2017) $
%    $Revision: 1233 $
%
% Copyright ?? 2011 The General Hospital Corporation (Boston, MA) "MGH"
%
% Terms and conditions for use, reproduction, distribution and contribution
% are found in the 'FreeSurfer Software License Agreement' contained
% in the file 'LICENSE' found in the FreeSurfer distribution, and here:
%
% https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
%
% Reporting: freesurfer@nmr.mgh.harvard.edu
%


  err = 1;

  if(nargin ~= 2 && nargin ~= 3)
    fprintf('USAGE: err = write_wfile(fname, w, <v>) \n');
    return;
  end

  vnum = length(w) ;

  % Handle when v not given or is empty %
  if (exist('v','var') ~= 1), v = []; end
  if (isempty(v)), v = 0:vnum-1; end

  % open it as a big-endian file
  fid = fopen(fname, 'wb', 'b') ;
  if(fid == -1)
    fprintf('ERROR: could not open %s\n',fname);
    return;
  end

  fwrite(fid, 0, 'int16') ;
  fwrite3(fid, vnum) ;
  for i=1:vnum
    fwrite3(fid, v(i)) ;          % vertex number (0-based)
    fwrite(fid,  w(i), 'float') ; % vertex value
  end

  fclose(fid) ;

  err = 0;
end

function [curv] = write_curv(fname, curv, fnum)
% [curv] = write_curv(fname, curv, fnum)
%
% writes a curvature vector into a binary file
%				fname - name of file to write to
%				curv  - vector of curvatures
%				fnum  - # of faces in surface.
%


%
% write_curv.m
%
% Original Author: Bruce Fischl
% CVS Revision Info:
%    $Author: gaser $
%    $Date: 2017-12-04 00:26:25 +0100 (Mon, 04 Dec 2017) $
%    $Revision: 1233 $
%
% Copyright (C) 2002-2007,
% The General Hospital Corporation (Boston, MA). 
% All rights reserved.
%
% Distribution, usage and copying of this software is covered under the
% terms found in the License Agreement file named 'COPYING' found in the
% FreeSurfer source code root directory, and duplicated here:
% https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
%
% General inquiries: freesurfer@nmr.mgh.harvard.edu
% Bug reports: analysis-bugs@nmr.mgh.harvard.edu
%

% assume fixed tetrahedral topology
if nargin == 2
  fnum = (length(curv)-2)*2;
end

% open it as a big-endian file
fid = fopen(fname, 'w', 'b') ;
vnum = length(curv) ;
NEW_VERSION_MAGIC_NUMBER = 16777215;
fwrite3(fid, NEW_VERSION_MAGIC_NUMBER ) ;
fwrite(fid, vnum,'int32') ;
fwrite(fid, fnum,'int32') ;
fwrite(fid, 1, 'int32');
fwrite(fid, curv, 'float') ;
fclose(fid) ;

end

function [curv, fnum] = read_curv(fname)
%
% [curv, fnum] = read_curv(fname)
% reads a binary curvature file into a vector
%
%
% read_curv.m
%
% Original Author: Bruce Fischl
% CVS Revision Info:
%    $Author: gaser $
%    $Date: 2017-12-04 00:26:25 +0100 (Mon, 04 Dec 2017) $
%    $Revision: 1233 $
%
% Copyright ?? 2011 The General Hospital Corporation (Boston, MA) "MGH"
%
% Terms and conditions for use, reproduction, distribution and contribution
% are found in the 'FreeSurfer Software License Agreement' contained
% in the file 'LICENSE' found in the FreeSurfer distribution, and here:
%
% https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
%
% Reporting: freesurfer@nmr.mgh.harvard.edu
%


%fid = fopen(fname, 'r') ;
%nvertices = fscanf(fid, '%d', 1);
%all = fscanf(fid, '%d %f %f %f %f\n', [5, nvertices]) ;
%curv = all(5, :)' ;

% open it as a big-endian file
if ~exist(fname,'file')
  error('cat_io_FreeSurfer:read_curv','Curvature file "%s" does not exist!\n', fname);
end
fid = fopen(fname, 'r', 'b') ;
if (fid < 0)
	 error('cat_io_FreeSurfer:read_curv','Could not open curvature file "%s"!\n', fname);
end
vnum = fread3(fid) ;
NEW_VERSION_MAGIC_NUMBER = 16777215;
if (vnum == NEW_VERSION_MAGIC_NUMBER)
	 vnum = fread(fid, 1, 'int32') ;
	 fnum = fread(fid, 1, 'int32') ;
   x    = fread(fid, 1, 'int32') ;
   curv = fread(fid, vnum, 'float') ; 
	   	
  fclose(fid) ;
else

	fnum = fread3(fid) ;
  curv = fread(fid, vnum, 'int32') ./ 100 ; 
  fclose(fid) ;
end

end

function write_annotation(filename, vertices, label, ct)
% Contact ythomas@csail.mit.edu or msabuncu@csail.mit.edu for bugs or questions 
%
%=========================================================================
%
%  Copyright (c) 2008 Thomas Yeo and Mert Sabuncu
%  All rights reserved.
%
%Redistribution and use in source and binary forms, with or without
%modification, are permitted provided that the following conditions are met:
%
%    * Redistributions of source code must retain the above copyright notice,
%      this list of conditions and the following disclaimer.
%
%    * Redistributions in binary form must reproduce the above copyright notice,
%      this list of conditions and the following disclaimer in the documentation
%      and/or other materials provided with the distribution.
%
%    * Neither the names of the copyright holders nor the names of future
%      contributors may be used to endorse or promote products derived from this
%      software without specific prior written permission.
%
%THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
%ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
%WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
%DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
%ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
%(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
%LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
%ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
%(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
%SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.    
%
%=========================================================================

  % write_annotation(filename, vertices, label, ct)
  %
  % Only writes version 2...
  %
  % vertices expected to be simply from 0 to number of vertices - 1;
  % label is the vector of annotation
  %
  % ct is a struct
  % ct.numEntries = number of Entries
  % ct.orig_tab = name of original ct
  % ct.struct_names = list of structure names (e.g. central sulcus and so on)
  % ct.table = n x 5 matrix. 1st column is r, 2nd column is g, 3rd column
  % is b, 4th column is flag, 5th column is resultant integer values
  % calculated from r + g*2^8 + b*2^16 + flag*2^24. flag expected to be all 0

  fp = fopen(filename, 'w', 'b');

  % first write vertices and label

  count = fwrite(fp, int32(length(label)), 'int');
  if(count~=1)
     error('write_annotation: Writing #vertices/labels not successful!!');
  end

  temp = zeros(length(label)*2,1);
  temp(1:2:end) = vertices;
  temp(2:2:end) = label;
  temp = int32(temp);

  count = fwrite(fp, int32(temp), 'int');
  if(count~=length(temp))
     error('write_annotation: Writing labels/vertices not successful!!');
  end

  %Write that ct exists
  count = fwrite(fp, int32(1), 'int');
  if(count~=1)
     error('write_annotation: Unable to write flag that ct exists!!');
  end

  %write version number
  count = fwrite(fp, int32(-2), 'int');
  if(count~=1)
      error('write_annotation: Unable to write version number!!');
  end

  %write number of entries
  count = fwrite(fp, int32(ct.numEntries), 'int');
  if(count~=1)
      error('write_annotation: Unable to write number of entries in ct!!');
  end

  %write original table
  orig_tab = [ct.orig_tab char(0)];
  count = fwrite(fp, int32(length(orig_tab)), 'int');
  if(count~=1)
      error('write_annotation: Unable to write length of ct source!!');
  end

  count = fwrite(fp, orig_tab, 'char');
  if(count~=length(orig_tab))
      error('write_annotation: Unable to write orig_tab!!');
  end

  %write number of entries
  count = fwrite(fp, int32(ct.numEntries), 'int');
  if(count~=1)
      error('write_annotation: Unable to write number of entries in ct!!');
  end

  %write ct
  for i = 1:ct.numEntries
      count = fwrite(fp, int32(i-1), 'int');
      if(count~=1)
          error('write_annotation: Unable to write structure number!!');
      end

      structure_name = [ct.struct_names{i} char(0)];
      count = fwrite(fp, int32(length(structure_name)), 'int');
      if(count~=1)
          error('write_annotation: Unable to write length of structure name!!');
      end
      count = fwrite(fp, structure_name, 'char');
      if(count~=length(structure_name))
          error('write_annotation: Unable to write structure name!!');
      end

      for j=1:4
          count = fwrite(fp, int32(ct.table(i, j)), 'int');
          if(count~=1)
             error('write_annotation: Unable to write red color'); 
          end
      end
      
  end

  fclose(fp);
end
function [vertices, label, colortable] = Read_Brain_Annotation(filename, varargin)
%
% NAME
%
%       function [vertices, label, colortable] = ...
%                                       read_annotation(filename [, verbosity])
%
% ARGUMENTS
% INPUT
%       filename        string          name of annotation file to read
%
% OPTIONAL
%       verbosity       int             if true (>0), disp running output
%                                       + if false (==0), be quiet and do not
%                                       + display any running output
%
% OUTPUT
%       vertices        vector          vector with values running from 0 to
%                                       + size(vertices)-1
%       label           vector          lookup of annotation values for 
%                                       + corresponding vertex index.
%       colortable      struct          structure of annotation data
%                                       + see below
%       
% DESCRIPTION
%
%       This function essentially reads in a FreeSurfer annotation file
%       <filename> and returns structures and vectors that together 
%       assign each index in the surface vector to one of several 
%       structure names.
%       
% COLORTABLE STRUCTURE
% 
%       Consists of the following fields:
%       o numEntries:   number of entries
%       o orig_tab:     filename of original colortable file
%       o struct_names: cell array of structure names
%       o table:        n x 5 matrix
%                       Columns 1,2,3 are RGB values for struct color
%                       Column 4 is a flag (usually 0)
%                       Column 5 is the structure ID, calculated from
%                       R + G*2^8 + B*2^16 + flag*2^24
%                       
% LABEL VECTOR
% 
%       Each component of the <label> vector has a structureID value. To
%       match the structureID value with a structure name, lookup the row
%       index of the structureID in the 5th column of the colortable.table
%       matrix. Use this index as an offset into the struct_names field
%       to match the structureID with a string name.      
%
% PRECONDITIONS
%
%       o <filename> must be a valid FreeSurfer annotation file.
%       
% POSTCONDITIONS
%
%       o <colortable> will be an empty struct if not embedded in a
%         FreeSurfer annotation file. 
%       

%
% read_annotation.m
% Original Author: Bruce Fischl
% CVS Revision Info:
%    $Author: gaser $
%    $Date: 2017-12-04 00:26:25 +0100 (Mon, 04 Dec 2017) $
%    $Revision: 1233 $
%
% Copyright ?? 2011 The General Hospital Corporation (Boston, MA) "MGH"
%
% Terms and conditions for use, reproduction, distribution and contribution
% are found in the 'FreeSurfer Software License Agreement' contained
% in the file 'LICENSE' found in the FreeSurfer distribution, and here:
%
% https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
%
% Reporting: freesurfer@nmr.mgh.harvard.edu
%

fp = fopen(filename, 'r', 'b');

verbosity = 0;
if ~isempty(varargin)
    verbosity       = varargin{1};  
end;

if(fp < 0)
   if verbosity, disp('Annotation file cannot be opened'); end;
   return;
end

A = fread(fp, 1, 'int');

tmp = fread(fp, 2*A, 'int');
vertices = tmp(1:2:end);
label = tmp(2:2:end);

bool = fread(fp, 1, 'int');
if(isempty(bool)) %means no colortable
   if verbosity, disp('No Colortable found.'); end;
   colortable = struct([]);
   fclose(fp);
   return; 
end

if(bool)
    
    %Read colortable
    numEntries = fread(fp, 1, 'int');

    if(numEntries > 0)
        
        if verbosity, disp('Reading from Original Version'); end;
        colortable.numEntries = numEntries;
        len = fread(fp, 1, 'int');
        colortable.orig_tab = fread(fp, len, '*char')';
        colortable.orig_tab = colortable.orig_tab(1:end-1);

        colortable.struct_names = cell(numEntries,1);
        colortable.table = zeros(numEntries,5);
        for i = 1:numEntries
            len = fread(fp, 1, 'int');
            colortable.struct_names{i} = fread(fp, len, '*char')';
            colortable.struct_names{i} = colortable.struct_names{i}(1:end-1);
            colortable.table(i,1) = fread(fp, 1, 'int');
            colortable.table(i,2) = fread(fp, 1, 'int');
            colortable.table(i,3) = fread(fp, 1, 'int');
            colortable.table(i,4) = fread(fp, 1, 'int');
            colortable.table(i,5) = colortable.table(i,1) + colortable.table(i,2)*2^8 + colortable.table(i,3)*2^16 + colortable.table(i,4)*2^24;
        end
        if verbosity
            disp(['colortable with ' num2str(colortable.numEntries) ' entries read (originally ' colortable.orig_tab ')']);
        end
    else
        version = -numEntries;
        if verbosity
          if(version~=2)    
            disp(['Error! Does not handle version ' num2str(version)]);
          else
            disp(['Reading from version ' num2str(version)]);
          end
        end
        numEntries = fread(fp, 1, 'int');
        colortable.numEntries = numEntries;
        len = fread(fp, 1, 'int');
        colortable.orig_tab = fread(fp, len, '*char')';
        colortable.orig_tab = colortable.orig_tab(1:end-1);
        
        colortable.struct_names = cell(numEntries,1);
        colortable.table = zeros(numEntries,5);
        
        numEntriesToRead = fread(fp, 1, 'int');
        for i = 1:numEntriesToRead
            structure = fread(fp, 1, 'int')+1;
            if (structure < 0)
              if verbosity, disp(['Error! Read entry, index ' num2str(structure)]); end;
            end
            if(~isempty(colortable.struct_names{structure}))
              if verbosity, disp(['Error! Duplicate Structure ' num2str(structure)]); end;
            end
            len = fread(fp, 1, 'int');
            colortable.struct_names{structure} = fread(fp, len, '*char')';
            colortable.struct_names{structure} = colortable.struct_names{structure}(1:end-1);
            colortable.table(structure,1) = fread(fp, 1, 'int');
            colortable.table(structure,2) = fread(fp, 1, 'int');
            colortable.table(structure,3) = fread(fp, 1, 'int');
            colortable.table(structure,4) = fread(fp, 1, 'int');
            colortable.table(structure,5) = colortable.table(structure,1) + colortable.table(structure,2)*2^8 + colortable.table(structure,3)*2^16 + colortable.table(structure,4)*2^24;       
        end
        if verbosity 
          disp(['colortable with ' num2str(colortable.numEntries) ' entries read (originally ' colortable.orig_tab ')']);
        end
    end    
else
    if verbosity
        disp('Error! Should not be expecting bool = 0');    
    end;
end

fclose(fp);


end
