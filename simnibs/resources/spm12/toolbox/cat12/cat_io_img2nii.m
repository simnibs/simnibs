function varargout=cat_io_img2nii(img,c,verb)
% ______________________________________________________________________
% Convert img/hdr images to nii.
%
%  varargout=cat_io_img2nii(img,c,verb)
%
%  img  .. list of *.img images 
%  c    .. delete original data (default=1)
%  verb .. display something (default=1)
% 
% ______________________________________________________________________
% $Revision: 764 $  $Date: 2015-11-17 14:11:53 +0100 (Di, 17 Nov 2015) $
  if ~exist('img','var') || isempty(img)
    img = spm_select(Inf  ,'img','Select img files'); 
  end
  if isa(img,'char'), img=cellstr(img); end
  if ~exist('c','var'), c=1; end
  if ~exist('verb','var'), verb=1; end

  n=numel(img);
  if verb, spm_progress_bar('Init',n,'Filtering','Volumes Complete'); end
  for i=1:n
    [pp,ff]=spm_fileparts(img{i}); niifile = fullfile(pp,sprintf('%s.nii',ff));
    if ~exist(niifile,'file')
      if verb, fprintf('%60s: ',img{i}); end; tic; 
      h = spm_vol(img{i});
      I = spm_read_vols(h);
      h.fname(end-2:end)='nii';
      spm_write_vol(h,I);
      if c==1, delete([img{i}(1:end-3) 'hdr']); delete([img{i}(1:end-3) 'img']); end
      if nargout==1, varargout{1}{i}=h.fname; end
    else
      if verb, fprintf('%60s: still exist! ',img{i}); end
      if nargout==1, varargout{1}{i}=niifile; end
    end
    if verb, spm_progress_bar('Set',i); fprintf('%5.2fs\n',toc); end
  end
  if verb, spm_progress_bar('Clear'); end
end