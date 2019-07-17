function varargout = cat_tst_qa_cleaner(data,opt)
%% _____________________________________________________________________
%  Estimate quality grades of given rating of one (or more) protocols
%  with 2 to 6 grads to separate passed, (unassignable) and failed 
%  images, by finding the first peak in the image quality histogram  
%  and using its width (standard deviation) in a limited range. 
%  If multiple protocols are used, than use the site variable opt.site 
%  and use the site depending output rths.
%
%  The passed range can be variated by opt.cf with lower values for harder 
%  and higher values for softer thresholds (more passed images), where 
%  opt.cf=1, describes a range that is similar to about 1% BWP noise that 
%  is equal to 5 rps.
%  ROC evaluation showed that opt.cf=0.72 allows the best separation of 
%  images without and with artifacts, but if the majority of your data 
%  include light artifacts (e.g. by movements in young children) that 
%  a softer weighing, e.g. opt.cf=2, is preferable (maximum is 4). 
%
%  Use the selftest with randomly generated data to get a first impression:
%    cat_tst_qa_cleaner('test')
% _____________________________________________________________________
%
%  This tool is still in development / undert test:
%   * the combination of different sites is not finished
%   * multiside output required a 'stacked' output
%
%  [Pth,rth,sq,rths,rthsc,sqs] = cat_tst_qa_remover(data[,opt])
%
%    Pth      .. global threshold for passed images 
%                (for odd grades this is in the middle of the unassignable)
%    rth      .. all global threshold(s) between grads
%    sq       .. estimated first peak and its std, where the std depend on
%                the number of grades!
%    rths     .. site depending thresholds between grads of each input 
%    rthsc    .. site depending thresholds between grads of each input 
%                (global corrected, removed further low quality data)
%    sqs      .. site depending first peaks and stds of passed data 
%
%    data     .. array of quality ratings or xml-files
%    opt      .. option structure
%     .grads  .. number of grads (2:6, default=6, see below)
%     .cf     .. factor for harder/softer thresholds (defaults=0.72)
%     .figure .. display histogramm with colored ranges of grads
%                 1 - use current figure
%                 2 - create new figure (default)
%                 3 - use one test figure (default in the selftest)
% _____________________________________________________________________
%
%  Grades:
%    2 grads:
%      P   passed
%      F   failed
%    3 grads:
%      P   passed
%      U   unassignable
%      F   failed
%    4 grads:
%      P+  clear passed 
%      P-  just passed
%      F+  just failed
%      F-  clear failed
%    5 grads:
%      P+  clear passed 
%      P-  just passed
%      U   unassignable
%      F+  just failed
%      F-  clear failed
%    6 grads (default):
%      P+  clear passed 
%      P   passed 
%      P-  just passed
%      F+  just failed
%      F   failed
%      F-  clear failed
%
% ______________________________________________________________________
% Robert Dahnke 
% Structural Brain Mapping Group
% University Jena
% ______________________________________________________________________
% $Id: cat_tst_qa_cleaner.m 1233 2017-12-03 23:26:25Z gaser $ 


  clear th; 
  if ~exist('opt','var'), opt = struct(); end
  def.cf        = 0.72;                 % normalization factor for rating 
  def.grads     = 6;                    % number of grads (default = 6)
  def.model     = 1;                    % model used for rating
  def.figure    = 2;                    % figure=2 for new/own figure
  def.smooth    = 0;                    % smoothing of output data
  def.siterf    = 1000000;              % round factor to identify similar resolution level 
  def.siteavgperc = [0.10 0.90];        % ?
  opt = cat_io_checkinopt(opt,def); 
  opt.cf = max( 0 , min( 4 , opt.cf )); % limit of cf
  
  % test options
  %opt.model = 2;
  %opt.grads = 6;
  
  % if no intput is given use SPM select to get some xml-files
  if ~exist('data','var') || isempty(data)
    data = cellstr(spm_select(inf,'XML','select qa XML-files',{},pwd,'^cat_.*')); 
  elseif ischar(data)
    data = cellstr(data);
  end
  if isempty(data) || (iscell(data) && all(cellfun('isempty',data)))
    if nargout>=1, varargout{1} = 3; end
    if nargout>=2, varargout{2} = 3; end
    if nargout>=3, varargout{3} = [2.5 0.5]; end
    if nargout>=4, varargout{4} = 3*ones(size(data)); end
    if nargout>=5, varargout{5} = 3*ones(size(data)); end
    if nargout>=6, varargout{6} = repmat([2.5 0.5],numel(data),1); end

    return;
  end
  if iscell(data) && ~strcmp(data{1},'test')
    fprintf('Load XML data');
    P = data; 
    xml = cat_io_xml(data,struct(),'read',1); clear data; 
    for di=1:numel(xml)
      opt.site(di,1) = xml(di).qualityratings.res_RMS; 
      data(di,1)     = xml(di).qualityratings.NCR; 
    end,
  end
  

  % --------------------------------------------------------------------
  % If a site variable is given (e.g. by the RMS resolution) then call
  % the cleanup for each subset. The threshold will be collected in a 
  % vector [markthss x opt.grads] with the same length as data. 
  % Nevertheless an average threshold will is estimated as average of 
  % the percentual range give by opt.siteavgperc with e.g. [0.1 0.9] to
  % concider 80% of the data.
  %  -------------------------------------------------------------------
  if isfield(opt,'site')
    if numel(opt.site)~=numel(data),
      error('cat_tst_qa_cleaner:numelsitedata','Numer of elements in data and opt.site have to be equal.\n');
    end
    opt.site = round(opt.site*opt.siterf)/opt.siterf; 
    sites    = unique(opt.site); 
    markth   = zeros(numel(sites),opt.grads-1); 
    markths  = zeros(numel(data),opt.grads-1); 
    siteth   = zeros(numel(data),2); 
    for si=1:numel(sites)
      sdatai = find(opt.site==sites(si));
      opts = opt; 
      opts = rmfield(opts,'site');
      opts.figure = 0; 
      [Sth,markth(si,:),out{1:4}] = cat_tst_qa_cleaner(data(sdatai),opts); %#ok<ASGLU>
      markths(sdatai,:) = repmat(markth(si,:),numel(sdatai),1); 
      siteth(sdatai,:)  = out{4}; 
    end
    % estimate global threshold
    markthss = sortrows(markth);
    th = cat_stat_nanmean(markthss(max(1,min(numel(sites),round(numel(sites)*opt.siteavgperc(1)))):...
                                   max(1,min(numel(sites),round(numel(sites)*opt.siteavgperc(2)))),:),1); 
    sd  = out{3}; 
    thx = out{4}; 
    % modify local rating based on the global one                 
    markths2 = markths;
    markths2 = min(markths2,1.2*repmat(th,size(markths2,1),1)); % higher thresholds even for sides with low rating 
    markths2 = max(markths2,0.8*repmat(th,size(markths2,1),1)); % lower  thresholds even for sides with high rating 
    d  = data; 

  else
    %  -----------------------------------------------------------------
    %  Simulate data, if no data is given by several normal distributed
    %  random numbers.
    %  -----------------------------------------------------------------
    if exist('data','var') && ~(iscell(data) && strcmp(data{1},'test'))
      d = data; 
      if numel(d)==0, 
        if nargout>=1, varargout{1} = nan; end
        if nargout>=2, varargout{2} = nan(1,opt.grads); end
        if nargout>=3, varargout{3} = nan(1,2); end
        if nargout>=4, varargout{4} = nan(size(data)); end
        if nargout>=5, varargout{5} = nan(size(data)); end
        if nargout>=6, varargout{6} = nan(size(data)); end
        return;
      end
    elseif iscell(data) && strcmp(data{1},'test')
      % Testcases with different quality ratings
      scans      = 100; % number of scans (per site) for simulation
      testcase   = round(rand(1)*10);
      randoffset = 0.5*randn(1,4);

      switch testcase
        case 0 % good quality, no outlier group
          d = [2.0 + randoffset(1) + 0.1*randn(1,round(scans*0.80)), ...
               2.5 + randoffset(2) + 0.3*randn(1,round(scans*0.15)), ...
               4.0 + randoffset(3) + 1.0*randn(1,round(scans*0.03)), ...
               5.0 + randoffset(4) + 1.0*randn(1,round(scans*0.02))];
         case 1 % good quality, with average outlier group
          d = [2.0 + randoffset(1) + 0.1*randn(1,round(scans*0.40)), ...
               2.5 + randoffset(2) + 0.3*randn(1,round(scans*0.40)), ...
               4.0 + randoffset(3) + 1.0*randn(1,round(scans*0.15)), ...
               5.0 + randoffset(4) + 1.0*randn(1,round(scans*0.05))];
        case 2 % good-average quality, with outlier group 
          d = [2.0 + randoffset(1) + 0.1*randn(1,round(scans*0.10)), ...
               2.5 + randoffset(2) + 0.3*randn(1,round(scans*0.50)), ...
               4.0 + randoffset(3) + 1.0*randn(1,round(scans*0.30)), ...
               5.0 + randoffset(4) + 1.0*randn(1,round(scans*0.10))];
        case 3 % good-average quality, without outlier group 
          d = [2.0 + randoffset(1) + 0.2*randn(1,round(scans*0.10)), ...
               2.5 + randoffset(2) + 0.3*randn(1,round(scans*0.50)), ...
               3.0 + randoffset(3) + 1.0*randn(1,round(scans*0.30)), ...
               4.0 + randoffset(4) + 1.0*randn(1,round(scans*0.10))]; 
        case 4 % average to low quality, with light falloff  
          d = [3.0 + randoffset(1) + 0.2*randn(1,round(scans*0.10)), ...
               3.5 + randoffset(2) + 0.3*randn(1,round(scans*0.50)), ...
               4.0 + randoffset(3) + 1.0*randn(1,round(scans*0.30)), ...
               5.0 + randoffset(4) + 1.0*randn(1,round(scans*0.10))];   
        case 5 % high to good quality, with light falloff  
          d = [1.0 + randoffset(1) + 0.2*randn(1,round(scans*0.10)), ...
               1.5 + randoffset(2) + 0.3*randn(1,round(scans*0.50)), ...
               2.0 + randoffset(3) + 1.0*randn(1,round(scans*0.30)), ...
               3.0 + randoffset(4) + 1.0*randn(1,round(scans*0.10))];   
        case 6 % high quality, no outlier
          d = [1.0 + randoffset(1) + 0.1*randn(1,round(scans*0.80)), ...
               1.5 + randoffset(2) + 0.3*randn(1,round(scans*0.13)), ...
               3.0 + randoffset(3) + 0.3*randn(1,round(scans*0.05)), ...
               5.0 + randoffset(4) + 0.3*randn(1,round(scans*0.02))];
        case 7 % good quality with second average peak 
          d = [2.0 + randoffset(1) + 0.1*randn(1,round(scans*0.30)), ...
               3.0 + randoffset(2) + 0.2*randn(1,round(scans*0.40)), ...
               4.0 + randoffset(3) + 0.5*randn(1,round(scans*0.10)), ...
               5.0 + randoffset(4) + 0.5*randn(1,round(scans*0.10))];   
        case 8 % good quality with second low quality peak 
          d = [1.0 + randoffset(1) + 0.1*randn(1,round(scans*0.50)), ...
               4.0 + randoffset(2) + 0.2*randn(1,round(scans*0.30)), ...
               4.0 + randoffset(3) + 0.5*randn(1,round(scans*0.10)), ...
               5.0 + randoffset(4) + 0.5*randn(1,round(scans*0.10))];    
        case 9 % good quality with second average and third low quality peak 
          d = [1.5 + randoffset(1) + 0.2*randn(1,round(scans*0.20)), ...
               3.0 + randoffset(2) + 0.3*randn(1,round(scans*0.20)), ...
               4.5 + randoffset(3) + 0.2*randn(1,round(scans*0.10)), ...
               2.0 + randoffset(4) + 0.8*randn(1,round(scans*0.50))];          
        case 10 % good quality with second average and third low quality peak 
          d = [1.5 + randoffset(1) + 0.1*randn(1,round(scans*0.10)), ...
               3.0 + randoffset(2) + 0.2*randn(1,round(scans*0.10)), ...
               4.5 + randoffset(3) + 0.2*randn(1,round(scans*0.10)), ...
               2.5 + randoffset(4) + 1.0*randn(1,round(scans*0.60))];           
      end

      % remove high quality outlier and set them to normal
      cor = max(1,median(d)-std(d)/2);
      md= d<(cor); d(md) = cor + 0.05*randn(1,sum(md));
      
      % set selftest figure
      opt.figure = 3; 
    end

    
  %  -------------------------------------------------------------------
  %  Models:
  %  I start with several ideas that all based on a similar idea: to 
  %  find the first peak that is given by the subset of images without
  %  inferences and to use the variance of this peak for further scaling
  %  of subsets for other grads. As far as IQR is already scaled, we 
  %  can limit the variance value ... e.g. the rating has an error of 
  %  0-2 rps (0.0-0.2 mark points) that is very low for high-quality data
  %  and higher for low-quality data. Due to our the general subdivion 
  %  of the rating scale in +,o, and - (e.g. B+,B,B-) we got a subrange 
  %  of 3.33 rps (1/3 mark points) that gives some kind of upper limit.
  %  -------------------------------------------------------------------
    thx = nan; sd = nan; th = zeros(1,opt.grads-1); 
    switch opt.model
      case 0
        % only global thresholding ... 
        % this is just to use the color bar output 
        thx = 3; 
        sd  = 1; 
        th  = 1.5:1:100;
        th(6:end) = []; 
      case 1 
        % kmeans model:
        % * estimate peaks based on the histogram
        % * mix the first and second peak until it fits to 30% of the data 
        %   or until the number of loops is similar the number of peaks 
        % * use the std give by one BWP noise level (0.5) to describe the 
        %   variance the passed interval.
        
        hx = hist(d,0.5:1:5.5);
        peaks = sum(hx>(max(hx)/5))*3;
        [thx,sdx] = kmeans3D(d,peaks); sdx = sdx./thx;
        for i=1:peaks
          if sum(d<thx(i))/numel(d) < 0.3
            thx(1) = cat_stat_nanmean(thx(1:2));
            sdx(1) = cat_stat_nanstd(d(d<thx(1)));
          end
        end
        sd    = 0.25 / (opt.grads/2) * opt.cf; % 0.5 = 1% BWP noise
        th(1) = thx(1) - sdx(1) + 2*sd(1); %- mean(sdx(1:min(3,numel(sdx)))) 
        for i = 2:opt.grads-1
          th(i) = th(i-1) + 2*sd(1); % 
        end
      case 2
        % similar to case 1, but with std optimization based on the data 
        % ... surprisingly the simple model 1 works better
        
        hx = hist(d,0.5:1:5.5); 
        %for i=1:1, hx(2:end-1) = cat_stat_nanmean(cat(1,hx(1:end-2),hx(2:end-1),hx(3:end)),1); end
        peaks = sum(hx>(max(hx)/5))*3;
        [thx,sdx] = kmeans3D(d,peaks); sdx = sdx./thx;
        for i=1:peaks
          %if numel(thx)>i && sum(d<thx(i))/numel(d) < 0.05
          %  thx(1) = []; sdx(1) = [];
          if sum(d<thx(i))/numel(d) < 0.3 %numel(thx)>i && 
            thx(1) = cat_stat_nanmean(thx(1:2)); 
            sdx(1) = cat_stat_nanstd(d(d<thx(1)));
          end
        end
        sdx(1) = cat_stat_nanstd(d(d<thx(1)));
        [thx,sdx] = kmeans3D(d(d<=(max([min(d),thx(1)+sdx(1)]))),3); thx=thx(2); sdx=sdx(2);  %sdx = sdx./thx;
        sd    = min(1/3,max(1/6,sdx(1))) / (opt.grads/2) * opt.cf; % 0.5 = 1% BWP noise*16
        th(1) = thx(1) - sdx(1) + 2*sd(1);
        for i = 2:opt.grads-1
          th(i) = th(i-1) + 2*sd(1); % 2*2/3*
        end
      
    end
    
    markths  = repmat(mean(th(floor(opt.grads/2):ceil(opt.grads/2))),size(data));
    markths2 = markths;
    siteth   = repmat([thx(1) sd],numel(data),1); 
  end
  
  
  
%%  ---------------------------------------------------------------------
%  Print:
%  This part is just for to plot a colorated histogram and the percents
%  of images in each group.
%  ---------------------------------------------------------------------
  if opt.figure
    if opt.figure==2
      f = figure;
      set(f,'color','w')
    elseif opt.figure==3
      f = findobj('type','figure','name','qa_cleaner_test');
      if isempty(f), figure('name','qa_cleaner_test'); else figure(f(1)); clf(f(1)); end
    end
    box on;
    
    %figure
    ss = 0.05; 
    [h,r]  = hist(d,0.5:ss:10.5); 
    for i=1:opt.smooth, h(2:end-1) = cat_stat_nanmean(cat(1,h(1:end-2),h(2:end-1),h(3:end)),1); end
    sh = 1; %sum(h);
    
    % background histogram (all data)
    %bar(r,h/sh,'facecolor',[0.8 0.8 0.8],'edgecolor','none');
    %fill(r,h/sh,[0.8 0.8 0.8],'edgecolor','none');
    hold on
    
    yl = [0 max(h)+1]; ylim(yl);
    % main grid
    for i=1.5:6,       plot([i i],ylim,'color',[0.8 0.8 0.8]); end
    switch numel(th)
      case 1
        hx = h; hx(r> th(1)+ss) = 0;        fill(r,hx/sh,[0.0  0.5  0.2],'edgecolor','none');  
        hx = h; hx(r<=th(1)-ss) = 0;        fill(r,hx/sh,[0.6  0.15 0.1],'edgecolor','none');  
        % main values 
        text(5,yl(2)*0.93,sprintf('%5.2f%% passed',sum(d< th(1))/numel(d)*100)          ,'color',[0.0  0.5  0.2]);
        text(5,yl(2)*0.85,sprintf('%5.2f%% failed',sum(d>=th(1))/numel(d)*100)          ,'color',[0.6  0.15 0.1]);
      case 2
        hx = h; hx(r>=th(1)+ss) = 0;        fill(r,hx/sh,[0.0  0.5  0.2],'edgecolor','none');  
        hx = h; hx(r<=th(1) | r>th(2)) = 0; fill(r,hx/sh,[0.85 0.75 0.3],'edgecolor','none');  
        hx = h; hx(r<=th(2)-ss) = 0;        fill(r,hx/sh,[0.6  0.15 0.1],'edgecolor','none');  
        % main values 
        text(5,yl(2)*0.93,sprintf('%5.2f%% passed',sum(d<th(1))/numel(d)*100)           ,'color',[0.0  0.5  0.2]);
        text(5,yl(2)*0.89,sprintf('%5.2f%% unassignable' ,sum(d>=th(1) & d<th(2))/numel(d)*100),'color',[0.85 0.75 0.3]);
        text(5,yl(2)*0.85,sprintf('%5.2f%% failed',sum(d>=th(2))/numel(d)*100)          ,'color',[0.6  0.15 0.1]);
      case 3
        % plot
        hx = h; hx(r>=th(1)+ss) = 0;           fill(r,hx/sh,[0.0  0.5  0.2],'edgecolor','none');  
        hx = h; hx(r<=th(1)-ss | r>th(2)) = 0; fill(r,hx/sh,[0.7  0.8  0.2],'edgecolor','none');  
        hx = h; hx(r<=th(2)-ss | r>th(3)) = 0; fill(r,hx/sh,[0.9  0.6  0.4],'edgecolor','none');  
        hx = h; hx(r<=th(3)-ss) = 0;           fill(r,hx/sh,[0.6  0.15 0.1],'edgecolor','none');  
        % main values
        text(5,yl(2)*0.93,sprintf('%5.2f%% passed',sum(d< th(2))/numel(d)*100),'color',[0   0.7  0]);
        text(5,yl(2)*0.88,sprintf('%5.2f%% failed',sum(d>=th(2))/numel(d)*100),'color',[0.8 0.0  0]);
        % detailed values
        text(5,yl(2)*0.75,sprintf('%5.2f%% passed+',sum(d< th(1))/numel(d)*100)          ,'color',[0.0  0.5  0.2]);
        text(5,yl(2)*0.70,sprintf('%5.2f%% passed-',sum(d>=th(1) & d<th(2))/numel(d)*100),'color',[0.7  0.8  0.2]);
        text(5,yl(2)*0.65,sprintf('%5.2f%% failed+',sum(d>=th(2) & d<th(3))/numel(d)*100),'color',[0.9  0.6  0.4]);
        text(5,yl(2)*0.60,sprintf('%5.2f%% failed-',sum(d>=th(3))/numel(d)*100)          ,'color',[0.6  0.15 0.1]);
      case 4
        % plot
        hx = h; hx(r>=th(1)+ss) = 0;           fill(r,hx/sh,[0.0  0.5  0.2],'edgecolor','none');  
        hx = h; hx(r<=th(1)-ss | r>th(2)) = 0; fill(r,hx/sh,[0.4  0.7  0.1],'edgecolor','none');  
        hx = h; hx(r<=th(2)-ss | r>th(3)) = 0; fill(r,hx/sh,[0.85 0.75 0.3],'edgecolor','none');  
        hx = h; hx(r<=th(3)-ss | r>th(4)) = 0; fill(r,hx/sh,[0.75 0.3  0.2],'edgecolor','none');  
        hx = h; hx(r<=th(4)-ss) = 0;           fill(r,hx/sh,[0.6  0.15 0.1],'edgecolor','none');  
        % main values 
        text(5,yl(2)*0.93,sprintf('%5.2f%% passed',sum(d<th(2))/numel(d)*100)           ,'color',[0.0  0.5  0.2]);
        text(5,yl(2)*0.89,sprintf('%5.2f%% check' ,sum(d>=th(2) & d<th(3))/numel(d)*100),'color',[0.85 0.75 0.3]);
        text(5,yl(2)*0.85,sprintf('%5.2f%% failed',sum(d>=th(3))/numel(d)*100)          ,'color',[0.7  0.0  0.0]);
        % detailed values
        text(5,yl(2)*0.75,sprintf('%5.2f%% passed+',sum(d<th(1))/numel(d)*100)           ,'color',[0.0  0.5  0.2]);
        text(5,yl(2)*0.71,sprintf('%5.2f%% passed-',sum(d>=th(1) & d<th(2))/numel(d)*100),'color',[0.4  0.7  0.1]);
        text(5,yl(2)*0.67,sprintf('%5.2f%% unassignable',sum(d>=th(2) & d<th(3))/numel(d)*100),'color',[0.85 0.75 0.3]);
        text(5,yl(2)*0.63,sprintf('%5.2f%% failed+',sum(d>=th(3) & d<th(4))/numel(d)*100),'color',[0.75 0.3  0.2]);
        text(5,yl(2)*0.59,sprintf('%5.2f%% failed-',sum(d>=th(4))/numel(d)*100)          ,'color',[0.6  0.15 0.1]);
      case 5
        % plot
        testbar=0; % it would be cool to use bars but they failed at least in MATLAB R2013 and killed the axis positions...
        if testbar==1
          hx = h; hx(r>=th(1)+ss) = 0;           bar(r,hx/sh,'facecolor',[0.0  0.5  0.2],'edgecolor','none','barwidth',1);
          hx = h; hx(r<=th(1)-ss | r>th(2)) = 0; bar(r,hx/sh,'facecolor',[0.4  0.7  0.1],'edgecolor','none','barwidth',1);  
          hx = h; hx(r<=th(2)-ss | r>th(3)) = 0; bar(r,hx/sh,'facecolor',[0.7  0.8  0.2],'edgecolor','none','barwidth',1);  
          hx = h; hx(r<=th(3)-ss | r>th(4)) = 0; bar(r,hx/sh,'facecolor',[0.9  0.6  0.4],'edgecolor','none','barwidth',1);  
          hx = h; hx(r<=th(4)-ss | r>th(5)) = 0; bar(r,hx/sh,'facecolor',[0.75 0.3  0.2],'edgecolor','none','barwidth',1);  
          hx = h; hx(r<=th(5)-ss) = 0;           bar(r,hx/sh,'facecolor',[0.6  0.15 0.1],'edgecolor','none','barwidth',1);      
        else
          hx = h; hx(r>=th(1)+ss) = 0;           fill(r,hx/sh,[0.0  0.5  0.2],'edgecolor','none');
          hx = h; hx(r<=th(1)-ss | r>th(2)) = 0; fill(r,hx/sh,[0.4  0.7  0.1],'edgecolor','none');  
          hx = h; hx(r<=th(2)-ss | r>th(3)) = 0; fill(r,hx/sh,[0.7  0.8  0.2],'edgecolor','none');  
          hx = h; hx(r<=th(3)-ss | r>th(4)) = 0; fill(r,hx/sh,[0.9  0.6  0.4],'edgecolor','none');  
          hx = h; hx(r<=th(4)-ss | r>th(5)) = 0; fill(r,hx/sh,[0.75 0.3  0.2],'edgecolor','none');  
          hx = h; hx(r<=th(5)-ss) = 0;           fill(r,hx/sh,[0.6  0.15 0.1],'edgecolor','none'); 
        end
        % main values 
        text(5,yl(2)*0.93,sprintf('%5.2f%% passed',sum(d<th(3))/numel(d)*100) ,'color',[0   0.7  0]);
        text(5,yl(2)*0.88,sprintf('%5.2f%% failed',sum(d>=th(3))/numel(d)*100),'color',[0.8 0.0  0]);
        % detailed values
        text(5,yl(2)*0.75,sprintf('%5.2f%% passed+',sum(d<th(1))/numel(d)*100)           ,'color',[0.0  0.5  0.2]);
        text(5,yl(2)*0.70,sprintf('%5.2f%% passed' ,sum(d>=th(1) & d<th(2))/numel(d)*100),'color',[0.4  0.7  0.1]);
        text(5,yl(2)*0.65,sprintf('%5.2f%% passed-',sum(d>=th(2) & d<th(3))/numel(d)*100),'color',[0.7  0.8  0.2]);
        text(5,yl(2)*0.60,sprintf('%5.2f%% failed+',sum(d>=th(3) & d<th(4))/numel(d)*100),'color',[0.9  0.6  0.4]);
        text(5,yl(2)*0.55,sprintf('%5.2f%% failed' ,sum(d>=th(4) & d<th(5))/numel(d)*100),'color',[0.75 0.3  0.2]);
        text(5,yl(2)*0.50,sprintf('%5.2f%% failed-',sum(d>=th(5))/numel(d)*100)          ,'color',[0.6  0.15 0.1]);
    end
    xlim([min(r),6.5]); 
    
    % subgrid
    for i=5/6:1/3:6.4, plot([i i],[0 0.03]*max(ylim),'color',[0.2 0.2 0.2]); end
        
    QMC   = cat_io_colormaps('marks+',17);
    color = @(QMC,m) QMC(max(1,min(size(QMC,1),round(((m-1)*3)+1))),:);
    
    
    % colored main grads
    FS = get(gca,'Fontsize')*1.3;
    set(gca,'XTick',0.5:1:6.5,'XTickLabel',{'100','90','80','70','60','50','40'},'TickLength',[0.02 0.02]);
    % further color axis objects...
    axA = copyobj(gca,gcf); axB = copyobj(axA,gcf); axC = copyobj(gca,gcf); 
    axD = copyobj(gca,gcf); axE = copyobj(gca,gcf); axF = copyobj(gca,gcf);
    % set colors...
    set(axA,'YTick',[],'XTickLabel',{},'XTick',1,'XColor',color(QMC,1),'Color','none','XTicklabel','A','TickLength',[0 0],'Fontsize',FS,'Fontweight','bold');
    set(axB,'YTick',[],'XTickLabel',{},'XTick',2,'XColor',color(QMC,2),'Color','none','XTicklabel','B','TickLength',[0 0],'Fontsize',FS,'Fontweight','bold');
    set(axC,'YTick',[],'XTickLabel',{},'XTick',3,'XColor',color(QMC,3),'Color','none','XTicklabel','C','TickLength',[0 0],'Fontsize',FS,'Fontweight','bold');
    set(axD,'YTick',[],'XTickLabel',{},'XTick',4,'XColor',color(QMC,4),'Color','none','XTicklabel','D','TickLength',[0 0],'Fontsize',FS,'Fontweight','bold');
    set(axE,'YTick',[],'XTickLabel',{},'XTick',5,'XColor',color(QMC,5),'Color','none','XTicklabel','E','TickLength',[0 0],'Fontsize',FS,'Fontweight','bold');
    set(axF,'YTick',[],'XTickLabel',{},'XTick',6,'XColor',color(QMC,6),'Color','none','XTicklabel','F','TickLength',[0 0],'Fontsize',FS,'Fontweight','bold');
    hold off; 
    
    if isfield(opt,'site') && numel(sites>1);
      title(sprintf('Histogram (cf=%0.2f) - global treshold for multisite output (n=%d)',opt.cf,numel(sites)),'Fontsize',FS);
    else
      title(sprintf('Histogram (cf=%0.2f)',opt.cf),'Fontsize',FS);
    end
    xlabel('IQR (rps)','Fontsize',FS); 
    ylabel('number of scans','Fontsize',FS); 
  end
  %%
  MarkColor = cat_io_colormaps('marks+',40); 
  if isfield(opt,'site') && numel(sites)>1, globcorr = ' (global corrected)'; else globcorr = ''; end
  if exist('P','var')
    files = P(data<=markths2(:,3)); 
    fprintf('PASSED%s: %0.2f%%\n',globcorr,numel(files)/numel(data)*100)
    if 0
      iqrs  = [xml(data<=markths2(:,3)).qualityratings];
      for fi=1:numel(files)
        cat_io_cprintf(MarkColor(max(1,round( iqrs(fi).IQR/9.5 * size(MarkColor,1))),:),'  %s\n',files{fi,1});
      end
    else
      
    end
    
    % bad files ...
    files = P(data>markths2(:,3) & data<=markths2(:,4)); 
    fprintf('FAILED+%s: %0.2f%%\n',globcorr,numel(files)/numel(data)*100)
    if 1
      iqrs  = [xml(data>markths2(:,3) & data<=markths2(:,4)).qualityratings];
      for fi=1:numel(files)
        cat_io_cprintf(MarkColor(max(1,round( iqrs(fi).IQR/9.5 * size(MarkColor,1))),:),'  %s\n',files{fi,1});
      end
    end
    files = P(data>markths2(:,4) & data<=markths2(:,5)); 
    iqrs  = [xml(data>markths2(:,4) & data<=markths2(:,5)).qualityratings];
    if 1
      fprintf('FAILED%s: %0.2f%%\n',globcorr,numel(files)/numel(data)*100)
      for fi=1:numel(files)
        cat_io_cprintf(MarkColor(max(1,round( iqrs(fi).IQR/9.5 * size(MarkColor,1))),:),'  %s\n',files{fi,1});
      end
    end
    files = P(data>markths2(:,5)); 
    fprintf('FAILED-%s: %0.2f%%\n',globcorr,numel(files)/numel(data)*100)
    if 1
      iqrs  = [xml(data>markths2(:,5)).qualityratings];
      for fi=1:numel(files)
        cat_io_cprintf(MarkColor(max(1,round( iqrs(fi).IQR/9.5 * size(MarkColor,1))),:),'  %s\n',files{fi,1});
      end
    end
  end
  
  
  %% create output
  if nargout>=1, varargout{1} = mean(th(floor(opt.grads/2):ceil(opt.grads/2))); end
  if nargout>=2, varargout{2} = th; end
  if nargout>=3, varargout{3} = [thx(1) sd(1)]; end
  if nargout>=4, varargout{4} = markths;  end
  if nargout>=5, varargout{5} = markths2; end
  if nargout>=6, varargout{6} = siteth; end
end
