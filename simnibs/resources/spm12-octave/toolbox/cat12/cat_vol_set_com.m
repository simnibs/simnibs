function cat_vol_set_com(vargin)
% use center-of-mass (COM) to roughly correct for differences in the
% position between image and template


if nargin == 1
	P = char(vargin.data);
else
  P = spm_select(Inf,'image','Select images to filter');
end
V = spm_vol(P);
n = size(P,1);

% pre-estimated COM of MNI template
com_reference = [0 -20 -15];

for i=1:n
  fprintf('Correct center-of-mass for %s\n',V(i).fname);
  Affine = eye(4);
  vol = spm_read_vols(V(i));
  avg = mean(vol(:));
  avg = mean(vol(vol>avg));
  
  % don't use background values
  [x,y,z] = ind2sub(size(vol),find(vol>avg));
  com = V(i).mat(1:3,:)*[mean(x) mean(y) mean(z) 1]';
  com = com';

  M = spm_get_space(V(i).fname);
  Affine(1:3,4) = (com - com_reference)';
  spm_get_space(V(i).fname,Affine\M);
end
