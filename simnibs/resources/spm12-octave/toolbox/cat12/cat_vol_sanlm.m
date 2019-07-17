function varargout = cat_vol_sanlm(varargin)
% Spatial Adaptive Non Local Means (SANLM) Denoising Filter
%_______________________________________________________________________
% Filter a set of images and add the prefix 'sanlm_'.
% Missing input (data) will call GUI or/and use defaults. 
%
% Input:
% job    - harvested job data structure (see matlabbatch help)
% 
% Output:
% out    - computation results, usually a struct variable.
%
% cat_vol_sanlm(job)
%   job.data    = set of images 
%   job.prefix  = prefix for filtered images (default = 'sanlm_') 
%   job.postfix = postfix for filtered images (default = '')
%   job.verb    = verbose processing (default = 1)
%   job.rician  = noise distribution
%   job.NCstr   = strength of noise correction (default = -inf) 
%     A value of 1 used the full correction by the SANLM filter. Values 
%     between 0 and 1 mix the original and the filtered images, whereas
%     INF estimated a global value depending on the changes of the SANLM
%     filter that reduce to strong filtering in high quality data.
%     Negative values work on a local rather than a global level. 
%     A value of -INF is recommend, but you define a set of values (see
%     example) for further changes depending on your data. Using multiple
%     values will add the string "[p|n]###_" to the prefix to describe the
%     ith input of job.NCstr, e.g. job.NCstr=[0.74,-0.33] will result in
%     "p074" and "n033".
%     
%     0                  .. no denoising 
%     1                  .. full denoising (original sanlm)
%     0 < job.NCstr < 1  .. automatic global correction with user weighting
%    -9 < job.NCstr < 0  .. automatic local correction with user weighting
%     inf                .. global automatic correction
%    -inf                .. local automatic correction
%
% Example:
%   cat_vol_sanlm(struct('data','','prefix','n','rician',0));
%   cat_vol_sanlm(struct('data','','NCstr',[-1:0.5:1,inf,-inf));
%_______________________________________________________________________
% Christian Gaser
% $Id: cat_vol_sanlm.m 1265 2018-02-01 10:10:55Z gaser $

if nargin == 0 
    job.data = cellstr(spm_select([1 Inf],'image','select images to filter'));
else
    job = varargin{1};
end

if nargout > 0
  varargout{1} = struct();
end

if ~isfield(job,'data') || isempty(job.data)
   job.data = cellstr(spm_select([1 Inf],'image','select images to filter'));
else
   job.data = cellstr(job.data);
end
if isempty(char(job.data)); return; end


def.verb    = 1;         % be verbose
def.prefix  = 'sanlm_';  % prefix
def.postfix = '';        % postfix
def.NCstr   = -Inf;      % 0 - no denoising, eps - light denoising, 1 - maximum denoising, inf = auto; 
def.rician  = 0;         % use inf for GUI
def.ares    = [1 1 2.5]; % resolution depending filter strength [active full-correction no-correction]
job = cat_io_checkinopt(job,def);


if strcmp(job.postfix,'NCstr'), job.postfix = sprintf('_NCstr%0.2f',job.NCstr); end
if ~isinf(job.NCstr), job.NCstr = max(-9.99,min(1,job.NCstr)); end           % guarantee values from -9.99 to 1 or inf
if isinf(job.rician), spm_input('Rician noise?',1,'yes|no',[1,0],2); end  % GUI


%%
V  = spm_vol(char(job.data));
Vo = V;
spm_clf('Interactive'); 
spm_progress_bar('Init',numel(job.data),'SANLM-Filtering','Volumes Complete');
for i = 1:numel(job.data)
    [pth,nm,xt,vr] = spm_fileparts(deblank(V(i).fname)); clear xt;  %#ok<ASGLU>

    stime = clock; 
    vx_vol  = sqrt(sum(V(i).mat(1:3,1:3).^2));

    src = single(spm_read_vols(V(i)));
    % prevent NaN and INF
    src(isnan(src) | isinf(src)) = 0;

    if numel(job.NCstr)>1
      fprintf('SANLM filtering: ');  
    end
    
    % histogram limit
    [hsrc,hval] = hist(src(:),10000); hp = cumsum(hsrc)./sum(hsrc); tol = 0.002;
    if min(src(:))~=0, src(src<hval(find(hp>tol,1,'first'))) = hval(find(hp>tol,1,'first')); end
    src(src>hval(find(hp<(1-tol),1,'last')))  = hval(find(hp<(1-tol),1,'last')); 
    
    % use intensity normalisation because cat_sanlm did not filter values below ~0.01 
    th = max(mean(src(src(:)>0)),abs(mean(src(src(:)<0))));
    src = (src / th) * 100; 
    cat_sanlm(src,3,1,job.rician);
    src = (src / 100) * th; 
    
    % measures changes
    srco   = single(spm_read_vols(V(i)));
    NCrate = cat_stat_nanmean(abs(src(:)/th - srco(:)/th)); 
    
    % set actual filter rate - limit later!
    % the factor 15 was estimated on the BWP 
    NCstr                         = job.NCstr; 
    NCstr(isinf(NCstr) & NCstr>0) = 15 * NCrate;   
    NCstr(isinf(NCstr) & NCstr<0) = -1;   
    stime2 = clock; 
    
   
    for NCstri = 1:numel(NCstr)
        stime3 = clock; 

        srco   = single(spm_read_vols(V(i)));
        
        if NCstr(NCstri)<0
        % adaptive local denoising 

            %% prepare local map
            % use less filtering for lowres data to avoid anatomical blurring ???
            NCs = max(eps,abs(src - srco)/th); 
            % preserve anatomical details by describing the average changes
            % and not the strongest - this reduce the ability of artefact
            % correction!
            stdNC = std(NCs(NCs(:)~=0)); 
            NCsm  = cat_vol_median3(NCs,NCs>stdNC,true(size(NCs))); % replace outlier
            [NCsr,resT2] = cat_vol_resize(NCsm,'reduceV',vx_vol,2,32,'meanm'); clear NCsm; 
            NCsr  = cat_vol_localstat(NCsr,true(size(NCs)),1,1);
            NCsr  = cat_vol_smooth3X(NCsr,2/mean(resT2.vx_volr));
            NCsr  = cat_vol_resize(NCsr,'dereduceV',resT2); 
            NCso  = NCs;  
            
            % no correction of local abnormal high values (anatomy)
            NCs = NCsr + (NCso>stdNC & NCso<=stdNC*4 & NCso>NCsr*2 & NCso<NCsr*16) .* (-NCsr); 
            NCs = cat_vol_smooth3X(NCs,2/mean(resT2.vx_vol));  
            % weighting
            NCs = NCs ./ mean(NCs(:)); NCs = max(0,min(1,NCs * (15*NCrate) * abs(NCstr(NCstri)))); 
            % volume dependency (lower filtering in lower resolution)
            if job.ares(1), NCs = NCs .* max(0,min(1,1-(mean(vx_vol)-job.ares(2))/diff(job.ares(2:3)))); end

            % heavy outlier / artefacts
            NCs = max(NCs,(NCso>stdNC*6) | (NCso>NCsr*12)); clear NCso; 
            
            
            %srco2 = srco.*(1-NCs) + src.*NCs;
            %% mix original and noise corrected image
            srco = srco.*(1-NCs) + src.*NCs; 
            NCstr(NCstri) = -cat_stat_nanmean(NCs(:)); 
            clear NCs;

        elseif NCstr(NCstri)==1
        % no adaption (original filter)
           srco   = src; 
        
        elseif NCstr(NCstri)>0
        % (adaptive) global denoising  

            %if job.ares(1), NCstr(NCstri) = NCstr(NCstri) .* max(0,min(1,1-(mean(vx_vol)-job.ares(2))/diff(job.ares(2:3)))); end
            NCstr(NCstri) = min(1,max(0,NCstr(NCstri)));
            
            % mix original and noise corrected image
            srco   = srco*(1-NCstr(NCstri)) + src*NCstr(NCstri); 
            
        else
        % no denoising ... nothing to do
            
        end
         
        
        if numel(job.NCstr)>1
          % add NCstr to the filename 
            if job.NCstr(NCstri)<0
                NCstring = sprintf('n%03d_',abs(job.NCstr(NCstri)*100)); 
            else
                NCstring = sprintf('p%03d_',abs(job.NCstr(NCstri)*100));
            end
            
            Vo(i).fname = fullfile(pth,[job.prefix NCstring nm job.postfix '.nii' vr]);            
            
            % display result and link
            if job.verb
                if NCstri==1
                  fprintf('%0.0fs\n',etime(stime2,stime));  
                end
                fprintf('NCstr = %- 5.2f > %4.2f,%5.0fs: Output %s\n',job.NCstr(NCstri),abs(NCstr(NCstri)),etime(clock,stime3),...
                  spm_file(Vo(i).fname,'link','spm_display(''%s'')'));
            end
        else
            Vo(i).fname = fullfile(pth,[job.prefix nm job.postfix '.nii' vr]);
            
            % display result and link
            if job.verb
                fprintf('NCstr = %- 5.2f > %4.2f,%5.0fs: Output %s\n',job.NCstr(NCstri),abs(NCstr(NCstri)),etime(clock,stime),...
                  spm_file(Vo(i).fname,'link','spm_orthviews(''Image'',''%s'')'));
            end
        end
        Vo(i).descrip = sprintf('%s SANLM filtered (NCstr=%-4.2f > %0.2f)',V(i).descrip,job.NCstr(NCstri),abs(NCstr(NCstri)));
       

        % use only float precision
        Vo(i).dt(1) = 16; 
        spm_write_vol(Vo(i), srco);
        spm_progress_bar('Set',i);
    end
end
spm_progress_bar('Clear');

