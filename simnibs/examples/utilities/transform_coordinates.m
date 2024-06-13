%%%
% Transform coordinates from MNI space to subject space and back
%
% Please add the SimNIBS "matlab" folder to your matlab path
% and run from the "ernie" directory of the example dataset
% visit www.simnibs.org for more information

% MNI position for M1 (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2034289/)
mni_M1 = [-37, -21, 58];

% Calculate the M1 position in the "ernie" model
ernie_M1 = mni2subject_coords(mni_M1, 'm2m_ernie/');
disp('Ernie M1 corrdinates:')
disp(ernie_M1)
% You can open the ernie_T1fs_conform.nii.gz and check if the coordinates match

% You can also do the inverse transformation (from subject to MNI)
mni_M1 = subject2mni_coords(ernie_M1, 'm2m_ernie/');
disp('Transforming coordinates back to MNI space:')
disp(mni_M1)
% Because SimNIBS uses non-linear transformations, you will not recover the
% original coordinates exactly
