function S=cat_vol_smooth3X(S,s,filter) 
% S=cat_vol_smooth3X(S,s,filter) 
% _________________________________________________________________________
% TODO - vx_vol!!!!!!!
%      - filter strength is also influcenced by the downsampling!!!
  if ~exist('s','var'), s=1; end
  if ~exist('filter','var'); filter='s'; end

  S(isnan(S(:)) | isinf(-S(:)) | isinf(S(:)))=0;                                          % correct bad cases
  
  SO=S;
  if size(size(S),2) == 2
    slice = 1;
    S = repmat(S,1,1,2*s+1);
  elseif size(size(S),2) == 1
    error('ERROR: cat_vol_smooth3X: Input S has to be a matrix or volume!'); 
  else
    slice = 0;
  end
  
  if s>0 && s<0.5
    S  = smooth3(S,'gaussian',3,0.5)*s + S*(1-s);
  elseif s>=0.5 && s<=1.0
    S  = smooth3(S,'gaussian',3,s);
  elseif s>1.0 && all(size(S)>6) 
    SR = reduceRes(S);           
    SR = cat_vol_smooth3X(SR,s/2); 
    S  = dereduceRes(SR,size(S)); 
  elseif s>=1.0 && any(size(S)<=6)
    S  = smooth3(S,'gaussian',5,s); 
  elseif s==0
    % nothing to do
  else
    S  = smooth3(S,'gaussian',1,s); 
%    error('ERROR: smooth3: s has to be greater 0'); 
  end
  switch filter
    case {'min'}, S = min(S,SO);
    case {'max'}, S = max(S,SO);
  end
  
  S(isnan(S(:)) | isinf(-S(:)) | isinf(S(:)))=0;    
  
  if slice 
    S = S(:,:,2); 
  end
end
function D=reduceRes(D,method)
  if ~exist('method','var'), method='linear'; end
  if mod(size(D,1),2)==1, D(end+1,:,:)=D(end,:,:); end
  if mod(size(D,2),2)==1, D(:,end+1,:)=D(:,end,:); end
  if mod(size(D,3),2)==1, D(:,:,end+1)=D(:,:,end); end
  
  [Rx,Ry,Rz] = meshgrid(single(1.5:2:size(D,2)),single(1.5:2:size(D,1)),single(1.5:2:size(D,3)));
  D=cat_vol_interp3f(imgExp(single(D),0),Rx,Ry,Rz,method);
end
function D=dereduceRes(D,sD,method)
  if ~exist('method','var'), method='linear'; end
  sD=sD/2+0.25;
  [Rx,Ry,Rz] = meshgrid(single(0.75:0.5:sD(2)),single(0.75:0.5:sD(1)),single(0.75:0.5:sD(3)));
  
  D=cat_vol_interp3f(single(imgExp(D,0)),Rx,Ry,Rz,method);
end
function D2=imgExp(D,d)
  if nargin<2, d=1; end
  if d>1, D=imgExp(D,d-1); end
  if d>0
    D2=zeros(size(D)+1,class(D));
    D2(1:end-1,1:end-1,1:end-1) = D; clear D; 
    for i=1:2
      D2(1:end,1:end,end) = D2(1:end,1:end,end-1);
      D2(1:end,end,1:end) = D2(1:end,end-1,1:end);
      D2(end,1:end,1:end) = D2(end-1,1:end,1:end);
    end
  else
    D2=D;
  end
end
function D2=imgDeExp(D,d)
  if nargin<2, d=1; end
  if d>0, D2=D(1:end-d,1:end-d,1:end-d); else D2=D; end
end