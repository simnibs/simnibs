function cat_vol_display_label(job)
% ______________________________________________________________________
% Skript for SPM display functions to show label overlays of one or 
% multiple atlas maps. 
% ______________________________________________________________________
% Robert Dahnke 
% Structural Brain Mapping Group
% University Jena
% ______________________________________________________________________
% $Id: cat_vol_display_label.m 846 2016-01-28 13:23:42Z gaser $ 

%
% based on cg_addtruecolorimage.m
%
% TODO:
%  - use xml colormaps if available (not defiend yet)
%  - display atlas rois of multiple atlas maps
%

  SVNid = '$Rev: 846 $';
  
  if ~exist('job','var'), job = struct(); end

  
  % check input files
  if ~isfield(job,'data')
    job.data = spm_select([2 13],'image','Select anatomical image and image to overlay');
  end
  job.data = cellstr(job.data);
  if numel(job.data)>13, job.data(13:end) = []; end

  
  % read volume
  V = spm_vol(char(job.data));
  
  
  % get caxis
  if ~isfield(job,'caxis')
    mni=nan(1,numel(V)); mxi=mni;
    for i=2:numel(V) 
      [mni(i),mxi(i)] = mn_mx_val(V(i));
    end
    mn = cat_stat_nanstat1d(mni,'nanmin');
    mx = cat_stat_nanstat1d(mxi,'nanmax');
    if mn < 0
      mn = min([-mn mx]);
      mx = mn;
      mn = -mn;
    end
    job.caxis = [mn;mx]; % spm_input('Minimum/Maximum','+1','e',[mn mx],2);
  end
 
  job.colormap = 'lines';
  if ~isfield(job,'colormap')
    job.colormapc = spm_input('Colormap',1,'e','jet');
  else
    if ischar(job.colormap)
      job.colormapc = eval(sprintf('%s(%d)',job.colormap,diff(job.caxis(:))));
    end
  end
  % add background in colormap (end of colormap!)
  if any(job.colormapc(end,1:3)~=0)
    job.colormapc = [job.colormapc;  0 0 0]; 
  end
  job.colormapc = max(0,min(1,job.colormapc + 0.2*(rand(size(job.colormapc))-0.5))); 
  mask = std(job.colormapc,1,2)<0.2; 
  job.colormapc(mask,:) = max(0,min(1,job.colormapc(mask,:) + 0.5*(rand(sum(mask==1),3)-0.5))); 
  
  
  % get transparency
  if ~isfield(job,'prop')
    job.prop = 0.2; % spm_input('Overlay transparecy (0..1)','+1','e',0.2,1);
  end
  
  
  %%
  spm_atlas('list','installed','-refresh');
  if size(job.data,1)<=2
    spm_image('init',V(1).fname); 
    spm_orthviews('AddContext',1);
  else
    spm_check_registration(repmat(V(1).fname,numel(V)-1,1));
  end
  
  fprintf(sprintf('%s',repmat('\b',1,73*2)));
  % remove SPM image/orthviews cmd-line link
  for i=2:numel(V)
    fprintf(sprintf('%s',repmat('\b',1,numel(sprintf('Display %s,1x',V(1).fname))))); % only the first anatomical image
  end
  
  % new banner
  spm('FnBanner',mfilename,SVNid);
  % all files
  if numel(V)>2
    files = sprintf('''%s''',job.data{1});
    for i=2:numel(V); files = sprintf('%s;''%s''',files,job.data{i}); end
    dispall = [' (' spm_file('all','link',...
        sprintf('cat_vol_display_label(struct(''data'',{{%s}},''colormap'',''%s''))', ...
        files,job.colormap)) ')  '];
  end
      
  global st; 
  % add colormaps
  for i=2:numel(V)
    % print new label based cmd-line link
    if i==2,     fprintf('Display ');                                  
    elseif i==3, fprintf('%s',dispall);
    else         fprintf('        '); 
    end
    
    fprintf('%s\n',spm_file(V(i).fname,'link',...
      sprintf('cat_vol_display_label(struct(''data'',{{''%s'';''%s''}},''colormap'',''%s''))', ...
      job.data{1},job.data{i},job.colormap)));

    if exist('mni','var')
      spm_orthviews('addtruecolourimage',i-1,V(i).fname,...
        job.colormapc,job.prop,mni(i),mxi(i));
    else
      spm_orthviews('addtruecolourimage',i-1,V(i).fname,...
        job.colormapc,job.prop,job.caxis(2),job.caxis(1));
    end
    vx   = sqrt(sum(V(i).mat(1:3,1:3).^2));
    spm_orthviews('resolution',min(vx(1:2)));
    spm_orthviews('redraw'); 
    
  end
  
  % only one label for all maps :/
  if size(job.data,1)<=2
    [pp,dsp] = fileparts(V(i).fname);
    dsp = spm_atlas('Load',['dartel_' dsp]);
    for ii=1:numel(st.vols)
      if ~isempty(st.vols{ii})
        st.vols{ii}.display = dsp;
      end
    end
    spm_ov_display('redraw');
  else
    for i=2:numel(V)
      [pp,dsp] = fileparts(V(i).fname);
      spm_orthviews('Caption',i-1,{dsp},'FontSize',12,'FontWeight','Bold');
    end
    
  end
end

%_______________________________________________________________________
function [mn,mx] = mn_mx_val(vol)
  mn = Inf;
  mx = -Inf;
  for i=1:vol.dim(3),
    tmp = spm_slice_vol(vol,spm_matrix([0 0 i]),vol.dim(1:2),0);
    imx = max(tmp(find(isfinite(tmp))));
    if ~isempty(imx),mx = max(mx,imx);end
    imn = min(tmp(find(isfinite(tmp))));
    if ~isempty(imn),mn = min(mn,imn);end
  end;
end