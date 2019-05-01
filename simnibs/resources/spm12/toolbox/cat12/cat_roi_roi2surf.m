function varargout = cat_roi_roi2surf(job)
% This function maps roi data to the surface, by setting the cdata of 
% all vertices of a ROI to the specified field.
%
%   cat_surf_roi2surf(job)
%
%   job.rdata .. csv - ROI files
%               xml - ROI files
% 
%   job.vars .. set of fieldnames
%
% ______________________________________________________________________
% Robert Dahnke
% $Id: cat_roi_roi2surf.m 1209 2017-11-07 16:20:27Z gaser $
  
  if nargin == 1
    def.verb = 1;
    def.usefsaverage  = 1; 
    def.assuregifti   = 1;
    def.fsaverage     = fullfile(spm('dir'),'toolbox','cat12','templates_surfaces','lh.central.freesurfer.gii');  
    def.inflated      = fullfile(spm('dir'),'toolbox','cat12','templates_surfaces','lh.inflated.freesurfer.gii');  
    def.dartelaverage = fullfile(spm('dir'),'toolbox','cat12','templates_surfaces','lh.central.Template_T1_IXI555_MNI152_GS.gii');    
    
    job = cat_io_checkinopt(job,def);
  else
    error('Only batch mode'); 
  end

  % volume ROIs
  % surface ROIs > ROI in Filename???

  switch job.surf
    case 'freesurfer', surf = job.fsaverage;
    case 'inflated',   surf = job.inflated;
    case 'dartel',     surf = job.dartelaverage;
    case 'subject',    surf = '';
  end
  
  % sides
  sides = {'lh.','rh.'}; 
  S = struct('vertices',[],'faces',[]); 
  if ~isempty(surf)
    [pps,ffs,ees] = spm_fileparts(surf);
    
    for si=1:2
      SX = gifti(fullfile(pps,[sides{si} ffs(4:end) ees])); 
      S(si).vertices = SX.vertices; 
      S(si).faces    = SX.faces;
      clear SX;
    end
  end
  
  spm_progress_bar('Init',numel(job.rdata{1}),...
    sprintf('ROI2Surface\n%s',numel(job.rdata{1})),'ROIs Completed'); 

  sname = cell(numel(job.rdata),1,1,numel(sides)); 
  for rfi=1:numel(job.rdata)
    [pp,ff,ee] = spm_fileparts(job.rdata{rfi});
    ffname  = textscan(ff,'%q','delimiter','_'); 

    if job.verb 
      fprintf('process "%s":\n',job.rdata{rfi});
    end
    
    %% first we need to load the RIO tables in a similar style 
    clear xml;
    switch ee
      case '.csv'
        % load csv data
        roiname = ffname{1}{2};
        subname = ffname{1}{3}; 
        xml.ROI.(roiname) = cat_io_csv(job.rdata{rfi}); 
      case '.xml'
        % read xml
        subname = ffname{1}{2}; 
        xml = cat_io_xml(job.rdata{rfi}); 
      otherwise 
        % error
    end
    
    % extract atlas names
    atlas = fieldnames(xml.ROI);

    
    %% load individual surfaces
    if strcmp(job.surf,'subject')
      % have to look for the resampled surface 
      % if it not exist we may resample the individual surface 
      % if this does not exist we print an error
      [pp2,ppl] = spm_fileparts(pp);
      surf = fullfile(pp2,strrep(ppl,'label','surf'),['lh.central.' subname '.gii']);

      [pps,ffs,ees] = spm_fileparts(surf);

      for si=1:2
        SX = gifti(fullfile(pps,[sides{si} ffs(4:end) ees])); 
        S(si).vertices = SX.vertices; clear SX;
      end
    end
    
    %%
    for ai = 1:numel(atlas)
      
      % load ROI data
      Proi = cat_vol_findfiles(fullfile(spm('dir'),'toolbox','cat12','atlases_surfaces'),['*' atlas{ai} '*']);
      for si=1:2
        [ppr,ffr,eer] = spm_fileparts(Proi{si});
        switch eer
          case '.annot'
             [vertices, S(si).rdata] = cat_io_FreeSurfer('read_annotation',Proi{si});
          case '.gii'
            SX = gifti(Proi{si}); 
            S(si).rdata = SX.cdata; clear SX;
          otherwise
            S(si).rdata = cat_io_FreeSurfer('read_surf_data',Proi{si});
        end
        S(si).rdata = uint16(S(si).rdata);
      end
      
      fields = xml.ROI.(atlas{ai})(1,3:end);
      switch job.fields
        case 'all'
          % nothing to do
        otherwise
          fields = setunion(fields,job.fields);

          % warning for unknown ROIs
      end

      %% mapping
      for si=1:numel(S)
        for fi = 1:numel(fields)
          %% 
          fid = find(cellfun('isempty',strfind(xml.ROI.(atlas{ai})(1,:),fields{fi}))==0);
          S(si).cdata = mapROI2surf(S(si).rdata,xml.ROI.(atlas{ai}),fid);

          % save data 
          [pp2,ppl] = spm_fileparts(pp);
          sname{rfi,ai,fi,si} = fullfile(pp2,strrep(ppl,'label','surf'),sprintf('%s%s-%s.ROI.%s.gii',sides{si},atlas{ai},fields{fi},subname));
          save(gifti(struct('vertices',S(si).vertices,'faces',S(si).faces,'cdata',S(si).cdata)),sname{rfi,ai,fi,si} );
          
          if job.verb
            fprintf('Output %s\n',spm_file(sname{rfi,ai,fi,si},'link','cat_surf_display(''%s'')'));
          end
        end        
      end

    end
    spm_progress_bar('Set',rfi);

  end
  spm_progress_bar('Clear');
  
  if nargout>0
    varargout{1} = sname; 
  end
end
function cdata = mapROI2surf(rdata,tab,fid)
  cdata = zeros(size(rdata),'single'); 
  for ri=2:size(tab,1)
    cdata(rdata==tab{ri,1}) = tab{ri,fid};
  end
end