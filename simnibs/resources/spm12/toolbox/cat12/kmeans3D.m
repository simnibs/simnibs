function [mu,su,nu]= kmeans3D(y,k)
% K-means clustering
% FORMAT mu = kmeans3D (y,k)
% 
% y          data 
% k          Number of components
%
% mu         vector of class means 
%
% modified version of
% spm_kmeans1.m 1143 2008-02-07 19:33:33Z spm $
%_______________________________________________________________________
% Christian Gaser
% $Id: kmeans3D.m 1127 2017-05-05 07:49:49Z gaser $

y=y(:)';
N=length(y);

% Spread seeds evenly according to CDF
[x,i]=sort(y);
seeds=[1,2*ones(1,k-1)]*N/(2*k);
seeds=ceil(cumsum(seeds));

last_i=ones(1,N);
mu=x(seeds);
su=zeros(size(mu));
nu=zeros(size(mu));

d = zeros(k,length(y));
for loops=1:1000,  
 for j=1:k,
   d(j,:)=(y-mu(j)).^2;
 end
 [tmp,i]=min(d);
 if sum(i-last_i)==0
   % If assignment is unchanged
   break;
 else
   % Recompute centres
   for j=1:k,
     mu(j)=mean(y(i==j));
     su(j)=std(y(i==j));
     nu(j)=sum(i==j)./numel(y(:));
   end
   last_i=i;
 end
end  
