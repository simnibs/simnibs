function FO = cat_io_handle_pre(F,pre,addpre,existfile)
% Remove all known cat prefix typs from a filename (and check if this file exist). 
  if nargin<4, existfile = 1; end
  if nargin<3, addpre    = 0; end

  [pp,ff,ee] = fileparts(F); 

  if ~addpre
    prefix{1} = {'r','m','e'};
    prefix{2} = {'ma','mb','mc','mv','mn','mi','ml', ...                     % modifications of T
                 'pa','pb','p0','p1','p2','p3','p4','p5','pf','pp' ...  % segmentations of T
                 'en','eb','ev','ei','ej','em', ...                     % error maps
                 't1','t2','t3','t4'};                                  % thickness and distaces
    prefix{3} = {'ra0','ra1','rw0','rw1','rwa','rwb','rwj', ...
                                   'rd0','rd1','rda','rdb','rdj', ...
                       'rcs','lcs','bcs','acs','hcs'};
    prefix{4} = {'ercs','elcs','ebcs','eacs','ehcs'};

    for pf=1:numel(prefix)
      if numel(ff)>pf+1 && any(strcmp(ff(1:pf),prefix{pf})) && ...
        (~existfile || exist(fullfile(pp,[ff(pf+1:end) ee]),'file'))
         FN = cat_io_handle_pre(fullfile(pp,[ff(pf+1:end) ee]),'',addpre,existfile); 
         if (~existfile || exist(FN,'file')), [ppn,ffn] = fileparts(FN); ff=ffn; end
      end
    end
  end
  
  FO = fullfile(pp,[pre ff ee]);
end