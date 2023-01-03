function [m, struct_names]=subject_atlas(m,  subdir, atlas_name)
% Transforms an atlas to subct space, and add it to the mesh m as a
% node_data field
%
% Parameters:
%   m: mesh
%   subdir: path to subject directory, example: 'm2m_ernie/'
%   atlas_name: 'a2009s', 'DK40' or 'HCP_MMP1'
%           'a2009s': Destrieux atlas (FreeSurfer v4.5, aparc.a2009s)
%           Cite: Destrieux, C. Fischl, B. Dale, A., Halgren, E. A sulcal 
%           depth-based anatomical parcellation of the cerebral cortex. 
%           Human Brain Mapping (HBM) Congress 2009, Poster #541
%
%           'DK40': Desikan-Killiany atlas (FreeSurfer, aparc.a2005s)
%           Cite: Desikan RS, S�gonne F, Fischl B, Quinn BT, Dickerson BC,
%           Blacker D, Buckner RL, Dale AM, Maguire RP, Hyman BT, Albert MS,
%           Killiany RJ. An automated labeling system for subdividing the
%           human cerebral cortex on MRI scans into gyral based regions of
%           interest. Neuroimage. 2006 Jul 1;31(3):968-80. 
%
%           'HCP_MMP1': Human Connectome Project (HCP) Multi-Modal Parcellation
%           Cite: Glasser MF, Coalson TS, Robinson EC, et al. A multi-modal 
%           parcellation of human cerebral cortex. Nature. 2016;536(7615):171-178.
% 
% Returns:
%   m: mesh, with annotation data added as node_data field
%   struct_names: names of the anatomical structures

% Guilherme B Saturnino, 2019
assert(exist(subdir, 'dir')==7, ['Could not find directory ' subdir])
assert(any(strcmp(atlas_name,{'a2009s', 'DK40', 'HCP_MMP1'})), ...
       ['invalid atlas name: ' atlas_name])

tmp_folder = tempname;
[status,result] = system([simnibs_cli_call('subject_atlas')...
                         ' -m ' subdir ' -a ' atlas_name ' -o ' tmp_folder]);
if status ~= 0
    rmdir(tmp_folder, 's');
    error('There was an error running subject_atlas:\n %s',result)
end
fn_lh = dir([tmp_folder filesep 'lh.*']);
[m, struct_names] = mesh_load_fsannot(m, [tmp_folder filesep fn_lh.name]);
rmdir(tmp_folder, 's');

end