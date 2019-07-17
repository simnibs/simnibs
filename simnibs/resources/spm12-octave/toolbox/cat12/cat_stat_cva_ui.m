function [CVA] = cat_stat_cva_ui(action,varargin)
% VOI extraction of adjusted data and CVA
% FORMAT [CVA] = cat_stat_cva_ui('specify',xSPM,SPM,CVA)
%
% xSPM   - structure containing specific SPM details
%     xSPM.Ic  - indice of contrast (in SPM.xCon)
% SPM    - structure containing generic analysis details
%
% CVA.contrast -  contrast name
% CVA.name     -  CVA name
% CVA.c        -  contrast weights
% CVA.X        -  contrast subspace
% CVA.Y        -  whitened and adjusted data
% CVA.X0       -  null space of contrast
%
% CVA.XYZ      -  locations of voxels (mm)
% CVA.xyz      -  seed voxel location (mm)
% CVA.VOX      -  dimension of voxels (mm)
%
% CVA.V        -  canonical vectors  (data)
% CVA.v        -  canonical variates (data)
% CVA.W        -  canonical vectors  (design)
% CVA.w        -  canonical variates (design)
% CVA.C        -  canonical contrast (design)
%
% CVA.chi      -  Chi-squared statistics testing D >= i
% CVA.df       -  d.f.
% CVA.p        -  p-values
%
% also saved in CVA_*.mat in the SPM working directory
%
% FORMAT [CVA] = cat_stat_cva_ui('results',CVA)
% Display the results of a CVA analysis
%__________________________________________________________________________
%
% This routine allows one to make inferences about effects that are
% distributed in a multivariate fashion or pattern over voxels. It uses
% conventional canonical variates (CVA) analysis (also know as canonical
% correlation analysis, ManCova and linear discriminant analysis).  CVA is
% a complement to MVB, in that the predictor variables remain the design
% matrix and the response variable is the imaging data in the usual way.
% However, the multivariate aspect of this model allows one to test for
% designed effects that are distributed over voxels and thereby increase
% the sensitivity of the analysis.
%
% Because there is only one test, there is no multiple comparison problem.
% The results are shown in term of the maximum intensity projection of the
% (positive) canonical image or vector and the canonical variates based on
% (maximally) correlated mixtures of the explanatory variables and data.
%
% CVA uses the generalised eigenvalue solution to the treatment and
% residual sum of squares and products of a general linear model. The
% eigenvalues (i.e., canonical values), after transformation, have a
% chi-squared distribution and allow one to test the null hypothesis that
% the mapping is D or more dimensional. This inference is shown as a bar
% plot of p-values.  The first p-value is formally identical to that
% obtained using Wilks' Lambda and tests for the significance of any
% mapping.
%
% This routine uses the current contrast to define the subspace of interest
% and treats the remaining design as uninteresting. Conventional results
% for the canonical values are used after the data (and design matrix) have
% been whitened; using the appropriate ReML estimate of non-sphericity.
%
% CVA can be used for decoding because the model employed by CVA does not
% care about the direction of the mapping (hence canonical correlation
% analysis). However, one cannot test for mappings between nonlinear
% mixtures of regional activity and some experimental variable (this is
% what the MVB was introduced for).
%
% References:
%
% Characterizing dynamic brain responses with fMRI: a multivariate
% approach. Friston KJ, Frith CD, Frackowiak RS, Turner R. NeuroImage. 1995
% Jun;2(2):166-72.
%
% A multivariate analysis of evoked responses in EEG and MEG data. Friston
% KJ, Stephan KM, Heather JD, Frith CD, Ioannides AA, Liu LC, Rugg MD,
% Vieth J, Keber H, Hunter K, Frackowiak RS. NeuroImage. 1996 Jun;
% 3(3):167-174.
%__________________________________________________________________________
% Copyright (C) 2008-2014 Wellcome Trust Centre for Neuroimaging

%
% $Id: cat_stat_cva_ui.m 1219 2017-11-19 23:28:59Z gaser $
%
% modified version of
% Karl Friston
% spm_cva_ui.m 6242 2014-10-14 11:16:02Z guillaume 


%-Get figure handles
%--------------------------------------------------------------------------
Finter = spm_figure('FindWin','Interactive');
try, spm_results_ui('Clear'); end
spm_input('!DeleteInputObj');

%-Review old analysis or proceed with a new one
%--------------------------------------------------------------------------
if ~nargin || isempty(action)
    action = spm_input('Canonical Variates Analysis','!+1','b', ...
        {'New Analysis','Results'}, char({'specify','results'}), 1);
    if strcmpi(action,'results'), varargin = {}; end
end

switch lower(action)
    
    case 'specify'
        %==================================================================
        %                      C V A  :  S P E C I F Y
        %==================================================================
        
        if nargin < 2     
          spmmatfile = spm_select(1,'^SPM\.mat$','Select SPM.mat');
          swd = spm_file(spmmatfile,'fpath');
          load(fullfile(swd,'SPM.mat'));
          SPM.swd = swd;
          cd(SPM.swd);

          % correct path for surface if analysis was made with different SPM installation
          if ~exist(SPM.xVol.G,'file')
            % check for 32k meshes
            if SPM.xY.VY(1).dim(1) == 32492 || SPM.xY.VY(1).dim(1) == 64984
              fsavgDir = fullfile(spm('dir'),'toolbox','cat12','templates_surfaces_32k');
            else
              fsavgDir = fullfile(spm('dir'),'toolbox','cat12','templates_surfaces');
            end
            [SPMpth,SPMname,SPMext] = spm_fileparts(SPM.xVol.G);
            SPM.xVol.G = fullfile(fsavgDir,[SPMname SPMext]);
          end

          %-Check the model has been estimated
          %--------------------------------------------------------------------------
          try
            XYZ  = SPM.xVol.XYZ;
          catch
    
            %-Check the model has been estimated
            %----------------------------------------------------------------------
            str = { 'This model has not been estimated.';...
              'Would you like to estimate it now?'};
            if spm_input(str,1,'bd','yes|no',[1,0],1)
              SPM = spm_spm(SPM);
            else
              SPM = []; xSPM = [];
              return
            end
          end
          
          [xSPM.Ic,xCon] = spm_conman(SPM,'T&F',Inf,...
               '    Select contrasts...',' ',1);

        else
          if nargin > 1, xSPM = varargin{1}; end
          if nargin > 2, SPM  = varargin{2}; end
          if nargin > 3, CVA  = varargin{3}; end
        end
        
        header = get(Finter,'Name');
        set(Finter,'Name','Canonical Variates analysis')
        
        %-Contrast specification
        %------------------------------------------------------------------
        con    = SPM.xCon(xSPM.Ic).name;
        c      = SPM.xCon(xSPM.Ic).c;
        c      = full(c);
        
        %-VOI specification
        %------------------------------------------------------------------
        name = con;
        name_img     = ['CVA_' strrep(name,' ','_')];
        name     = ['CVA_' strrep(name,' ','_') '.mat'];
        
        xyzmm = [0 0 0];

        %-Specify search volume
        %------------------------------------------------------------------
        try
            xY = CVA.xY;
            CVA = rmfield(CVA,'xY');
        catch
            xY = [];
        end
        xY.xyz = xyzmm;
        
        Q      = ones(1,size(SPM.xVol.XYZ, 2));
  
        XYZmm  = SPM.xVol.M(1:3,:)*[SPM.xVol.XYZ; Q];
        XYZ = XYZmm;
        j = 1:size(XYZmm,2);
        
        %-Extract required data from results files
        %==================================================================
        
        spm('Pointer','Watch')
        
        %-Get explanatory variables (data)
        %------------------------------------------------------------------
        Y    = spm_data_read(SPM.xY.VY,'xyz',SPM.xVol.XYZ(:,j));
        
        if isempty(Y)
            spm('alert*',{'No voxels in this VOI';'Please use a larger volume'},...
                'Canonical Variates analysis');
            return
        end
        
        %-Remove serial correlations and get design (note X := W*X)
        %------------------------------------------------------------------
        Y   = SPM.xX.W*Y;
        X   = SPM.xX.xKXs.X;
        
        %-Null-space
        %------------------------------------------------------------------
        X0  = [];
        try, X0 = [X0 blkdiag(SPM.xX.K.X0)]; end          %-drift terms
        try, X0 = [X0 spm_detrend(SPM.xGX.gSF)]; end      %-global estimate
        
        
        %-Canonical Variate Analysis
        %==================================================================
%        U     = spm_mvb_U(Y,'compact',spm_svd([X0, X-X*c*pinv(c)]),XYZmm);
%        CVA   = spm_cva(Y, X, X0, c, U);
        CVA   = spm_cva(Y, X, X0, c);
        
        % check for inverted effects
        if any(sign(CVA.C)-sign(c))
          CVA.C = -CVA.C;
          CVA.V = -CVA.V;
        end
        
        % scale canonical variates to SD of 1
        for j=1:size(CVA.V,2)
          v = CVA.V(:,j);
          CVA.V(:,j) = v./(eps + std(v(v ~= 0 & ~isnan(v))));;
        end
        
        %-Save results
        %==================================================================
        M     = SPM.xVol.M(1:3,1:3);             %-voxels to mm matrix
        VOX   = sqrt(diag(M'*M))';               %-voxel dimensions
        
        %-Assemble results
        %------------------------------------------------------------------
        CVA.contrast = con;                      %-contrast name
        CVA.name     = name;                     %-CVA name
        
        CVA.XYZ      = XYZmm;                    %-locations of voxels (mm)
        CVA.xyz      = xyzmm;                    %-seed voxel location (mm)
        CVA.VOX      = VOX;                      %-dimension of voxels (mm)
        if spm_mesh_detect(SPM.xY.VY)
          CVA.xVol   = SPM.xVol;
        end
        
        %-Save
        %------------------------------------------------------------------
        save(fullfile(SPM.swd,name),'CVA', spm_get_defaults('mat.format'));
        
        for j=1:size(CVA.V,2)
          if CVA.p(j) < 0.05
            if spm_mesh_detect(SPM.xY.VY)
              F = spm_file(sprintf('%s_%02d',name_img,j),'ext','.gii');
              F     = spm_file(F,'CPath');
              fprintf('Canonical image %s will be saved (significance: P<%g)\n',F,CVA.p(j));
              M     = gifti(SPM.xVol.G);
              C     = zeros(1,size(SPM.xVol.G.vertices,1));
              C(SPM.xVol.XYZ(1,:)) = CVA.V(:,j); % or use NODE_INDEX
              M.cdata = C;
              save(M,F);
              cmd   = 'spm_mesh_render(''Disp'',''%s'')';
            else

              Y      = zeros(SPM.xVol.DIM(1:3)');
              OFF    = SPM.xVol.XYZ(1,:) + SPM.xVol.DIM(1)*(SPM.xVol.XYZ(2,:)-1 + SPM.xVol.DIM(2)*(SPM.xVol.XYZ(3,:)-1));
              Y(OFF) = CVA.V(:,j);

              VO = SPM.xY.VY(1);
              VO.fname = spm_file(sprintf('%s_%02d',name_img,j),'ext','.nii');
              fprintf('Canonical image %s will be saved (significance: P<%g)\n',VO.fname,CVA.p(j));
              VO.dt = [spm_type('float32') spm_platform('bigend')];
              spm_write_vol(VO,Y);

              cmd = 'spm_image(''display'',''%s'')';
            end
          end
        end
          
        assignin('base','CVA',CVA);
        
        %-Display results
        %------------------------------------------------------------------
        cat_stat_cva_ui('results',CVA);
        
        %-Reset title
        %------------------------------------------------------------------
        set(Finter,'Name',header)
        spm('Pointer','Arrow')
        
        
    case 'results'
        %==================================================================
        %                      C V A  :  R E S U L T S
        %==================================================================
        
        %-Get CVA if necessary
        %------------------------------------------------------------------
        if isempty(varargin)
            [CVA,sts] = spm_select(1,'mat',...
                'Select CVA to display',[],[],'^CVA.*\.mat$');
            if ~sts, return; end
        else
            CVA = varargin{1};
        end
        if ischar(CVA)
            CVA = load(CVA);
            CVA  = CVA.CVA;
        end
        
        %-Show results
        %------------------------------------------------------------------
        Fgraph = spm_figure('GetWin','Graphics');
        
        %-Unpack
        %------------------------------------------------------------------
        VOX      = CVA.VOX;
        XYZ      = CVA.XYZ;
                
        %-Maximum intensity projection (first canonical image)
        %------------------------------------------------------------------
        if isfield(CVA,'xVol');
          hMIPax = axes('Parent',Fgraph,'Position',[0.05 0.60 0.45 0.26],'Visible','off');
          hMax = cat_surf_render('Disp',CVA.xVol.G,'Parent',hMIPax);
          tmp = zeros(1,prod(CVA.xVol.DIM));
          tmp(CVA.xVol.XYZ(1,:)) = CVA.V(:,1);
          mx = max(abs(min(CVA.V(:,1))),abs(max(CVA.V(:,1))));
          hMax = cat_surf_render('Overlay',hMax,tmp);
          hMax = cat_surf_render('Colourbar',hMax,'on');
          hMax = cat_surf_render('Colourmap',hMax,jet(64));
          hMax = cat_surf_render('Clip',hMax,[-2 2]);
          hMax = cat_surf_render('Clim',hMax,[-mx mx]);
        else
          subplot(2,2,1)
          spm_mip(CVA.V(:,1).*(CVA.V(:,1) > 0),XYZ(1:3,:),diag(VOX));
          axis image
        end
        title({'(Principal) canonical image',[CVA.name ':' CVA.contrast]})

        %-Inference and canonical variates
        %------------------------------------------------------------------
        Xstr{1} = 'Dimensionality';
        Xstr{2} = ['Chi-squared: ' sprintf('%6.1f ',  CVA.chi)];
        Xstr{3} = ['           df: ' sprintf('%6.0f ',CVA.df) ];
        
        subplot(2,2,2)
        bar(log(CVA.p)); hold on
        plot([0 (length(CVA.p) + 1)],log(0.05)*[1 1],'r:','LineWidth',4), hold off
        xlabel(Xstr)
        ylabel('log p-value')
        axis square
        title({'Test of dimensionality';sprintf('minimum p = %.2e',min(CVA.p))})
        
        subplot(2,2,3)
        plot(CVA.w,CVA.v,'.')
        xlabel('prediction')
        ylabel('response')
        axis square
        title('Canonical variates')
        
        %-Canonical contrast
        %------------------------------------------------------------------
        i       = find(CVA.p < 0.05);
        str     = 'Significant canonical contrasts';
        if isempty(i)
            i   = 1;
            str = 'first canonical contrast';
        end
        subplot(2,2,4)
        bar(CVA.C(:,i))
        xlabel('Parameter')
        axis square
        title(str)        
        
    otherwise
        error('Unknown action.');
        
end
