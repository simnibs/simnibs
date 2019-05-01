function I=cat_vol_nanmean3(I,s,iterations)
% _________________________________________________________________________
% smooth image with nans
% _________________________________________________________________________
  if ~isa('I','single'), I = single(I); end; 
  if ~exist('s','var'), s=1; end
  if ~exist('iterations','var'), iterations=1; end
  for iteration=1:iterations
    I2 = I; I3 = I;
    for i=1+s:size(I,1)-s, I2(i,:,:) = cat_stat_nanmean(I3(i-s:i+s,:,:),1); end
    for i=1+s:size(I,2)-s, I3(:,i,:) = cat_stat_nanmean(I2(:,i-s:i+s,:),2); end
    for i=1+s:size(I,3)-s, I2(:,:,i) = cat_stat_nanmean(I3(:,:,i-s:i+s),3); end  
    I(isnan(I)) = I2(isnan(I));     
  end
end 