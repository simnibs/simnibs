function cat_tst_staple_multilabels(P,Pm,Q,verb)

  if nargin == 0
    P = spm_select(Inf,'image','Select images from different experts');
    Pm = spm_select(1,'image','Select mask image if needed');
    [pp,ff,ee] = spm_fileparts(P(1,:));
    Q = fullfile(pp,[ff '_gt' ee]);
  end
  if ~exist('verb','var'), verb=0; end

  n = size(P,1);

  V  = spm_vol(deblank(P));
  Vm = spm_vol(Pm);
  Vo = V(1);
  Vo.fname = Q;

  vol = zeros([V(1).dim(1:3) n],'single');
  for i=1:n
    vol(:,:,:,i) = single(spm_read_vols(V(i)));
  end 
  n=n+1; vol(:,:,:,n) = cat_stat_nanmean(vol,4);
  %vol(isinf(vol(:)) | isnan(vol(:))) = 0;
  vol_size = size(vol);

  % check dimensions for data+mask
  if ~isempty(Vm); V(n+1) = Vm; end
  if any(any(diff(cat(1,V.dim),1,1),1))
    error('Dimensions differ');
  end 

  % mask images and remove background
  if ~isempty(Vm)
    mask = spm_read_vols(Vm);
  else
    mask = zeros(V(1).dim(1:3),'single');
    for i=1:n
      volt = single(vol(:,:,:,i)>0.5);
      if sum(volt(:))>10000 % to avaid using of missclassified images for masking
        mask = mask + volt;
      end
    end
    mask = smooth3(single(mask))>=(max(mask(:))/2);
  end
  mask_ind = find(mask > 0);
  masked_vol = zeros(length(mask_ind),n,'single');
  for i=1:n
    tmp = vol(:,:,:,i); 
    tmp(isnan(tmp(:)) | isinf(tmp(:)) | tmp(:)<0)=0;
    masked_vol(:,i) = round(tmp(mask_ind));
  end
  meanvol=mean(vol,4);
  clear tmp vol;
  
  numLabel = max(masked_vol(:)) + 1;
  minLabel = min(masked_vol(:));
  if minLabel ~= 0
    error('Minimum value in data is not zero');
  end

  [W, Theta] = STAPLE_multiLabels_nD(masked_vol, numLabel,verb);

  slabel = zeros(size(mask),'single');
  slabel(mask_ind) = W - 1;
  slabel(slabel>meanvol*2)=0;
  Vo.dt(1) = 2; 
  spm_write_vol(Vo,slabel);
end

% Program:  function [W,Theta,stop]=STAPLE_multiLabels(expertSegmentations,numLabels)
% Date:     2005/01/25 
% Language: Matlab
% 
% AUTHORS:  Meritxell Bach Cuadra (http://ltswww.epfl.ch/~bach)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% REFERENCE FOR THE IMPLEMENTATION:
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% "Simultaneous Truth and Performance Level Estimation (STAPLE): An
% algorithm for the Validation of Image Segmentation", Warfield et al, 
% IEEE Transactions on Medical Imaging, Volume: 23 , Issue: 7 , July 2004 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TERMS OF USE: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% You can use/modify this program for any use you wish, provided you cite 
% the above references in any publication about it. 
%     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DISCLAIMER:
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  In no event shall the authors or distributors be liable to any party for 
%  direct, indirect, special, incidental, or consequential damages arising out 
%  of the use of this software, its documentation, or any derivatives thereof, 
%  even if the authors have been advised of the possibility of such damage.
%    
%  The authors and distributors specifically disclaim any warranties, including,
%  but not limited to, the implied warranties of merchantability, fitness for a
%  particular purpose, and non-infringement.  this software is provided on an 
%  "as is" basis, and the authors and distributors have no obligation to provide
%  maintenance, support, updates, enhancements, or modifications.
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to calculate the new label map
% using a EM-HMRF assuming Gaussian distribution.
%
% It will first use the normal EM algorithm to find an estimation
% of the means and sigmas of the tissue type:
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Usage:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dij=expertSegmentations: numLine x numCol x numExperts matrix (for the moment only 2D segmentations) of
% integer values (typically 0,1,2,...,L, where L is the number of different labels)
% 
% numLabels = number of present labels including the background
% 
% W is the estimated True Segmentation probability, it is a matrix of N x
% numLabels elements (where N is the total number of voxels in the image)
%
% Theta is a numLabels x numLabels x numExperts matrix. For each expert, a
% confusion table corresponding to the probability of an expert to classify
% the voxel as s when is actually s_prime (float values)
% 
% stop is a vector having the convergence criteria history (here the normalized 
% trace of the expert parameters is considered)

function [W,Theta,stop]=STAPLE_multiLabels_nD(expertSegmentations,numLabels,verb)

  % expertSegmentations will be converted to a 2D matrix: N x numExperts
  origsize = size(expertSegmentations);
  ndim = length(origsize) - 1;
  expertSegmentations = 1 + reshape(expertSegmentations,prod(origsize(1:ndim)),origsize(ndim+1));

  s=size(expertSegmentations);

  N=s(1);
  numExperts=s(2);

  % Prior probability (f(Ti=s)) of belonging to one class (sample mean of the relative proportion of the label in the segmentations)
  % Gamma_s=f(Ti=s)
  Gamma_s=zeros(numLabels,1,'single'); 
  for k=1:numLabels, 
      temp=zeros(s);
      temp(find(expertSegmentations==(k-1)))=1;
      Gamma_s(k)=sum(sum(sum(temp)))./(N*numExperts);
  end
  if verb
    fprintf('Prior probability (f(Ti=s)) of belonging to one class (Gamma):\n');
    fprintf('%5.5f  ',Gamma_s);
    fprintf('\n\n');
  end
  
  clear temp

  %Expert parameters: confusion table Theta
  Theta=zeros(numLabels,numLabels,numExperts,'single');
  % Initialization as proposed by Warfield: high value equal for all experts
  for j=1:numExperts,
      for s=1:numLabels,
          for k=1:numLabels,
              if(s==k),
                  Theta(s,k,j)=0.99;  
              else
                  Theta(s,k,j)=(1-0.99)./(numLabels-1);
              end
          end
      end
  end

  %Convergence
  epsilon=10^(-7); % As proposed in the paper
  maxIter=25;
  stop=zeros(maxIter,1);
  %Loop until convergence or maximum number of iterations
  k=1;
  while(k<=maxIter),
      if verb, fprintf('Iteration %d\n',k); end
      %Estimation of W
      W=zeros(N,numLabels,'single');
      %Eq.20 TMI
      sumW=zeros(N,1,'single');
      for s=1:numLabels,
          %numeratorW=zeros(N,1);
          numeratorW=Gamma_s(s);
          for j=1:numExperts,
              numeratorW=numeratorW.*Theta(expertSegmentations(:,j),s,j);
          end
          sumW=sumW+numeratorW;
          W(:,s)=numeratorW;
      end
      %Normalization among all labels
      for s=1:numLabels,
          W(:,s)=W(:,s)./sumW(:);
          %W(find(sumW==0),s)=0;
      end

      %Estimation of Theta (expert parameters)
      for j=1:numExperts,
          for s=1:numLabels,
              for s_prime=1:numLabels,
                  index=find(expertSegmentations(:,j)==s_prime);
                  temp=zeros(N,1);
                  temp(index)=W(index,s);
                  if(sum(W(:,s))~=0)
                      Theta(s_prime,s,j)=sum(temp)./sum(W(:,s));
                  end
              end
          end
      end
    clear temp

      %Compute stopping criteria
      for j=1:numExperts,
          stop(k)=(stop(k)+trace(Theta(:,:,j)));
      end
      stop(k)=stop(k)./numLabels/numExperts;
      if(k>1),
          if((stop(k-1)-stop(k))<epsilon),
              convergenceAt=k;
              k=maxIter;
          end
      end
      k=k+1;
  end

  Wnew = zeros(N,1,'single');
  for s=1:numLabels
    Wnew = Wnew + single(s)*W(:,s);							   
  end

  if ndim > 1
    W = reshape(Wnew,origsize(1:ndim));
  else
    W = Wnew;
  end

  % Showing convergence history
  if verb
    figure;
    index=find(stop);
    plot(stop(index))
    title('STAPLE convergence evolution')
    xlabel('Iteration')
    ylabel('Normalized trace of expert parameters')
  end
end