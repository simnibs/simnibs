%
% example script that reads simulation results mapped onto the
% freesurfer surfaces of the subject. It extracts peak magnitude and some
% percentiles of E. This is done (i) for the whole cortex, (ii) separately
% for left and right hemisphere, and (iii) in all regions defined in the
% aparc.a2009s annotation files
%
% A. Thielscher, 06-Nov-2017


% input files and options
%-------------------------

% get results for both hemispheres, or only for the one specified
% by the surface name below
bothHemi=true;

% name of a freesurfer or CAT12 surface
%
% freesurfer surfaces: fs_{subID}/surf/lh.pial, fs_{subID}/surf/rh.pial
% CAT12 surfaces: m2m_{subID}/segment/cat/surf/lh.central.T1fs_conform.gii, m2m_{subID}/segment/cat/surf/rh.central.T1fs_conform.gii
fname_surface='.../m2m_test/segment/cat/surf/lh.central.T1fs_conform.gii'; % UPDATE

% corresponding name of a results file in the 'subject_overlays' subfolder
%
% NOTE:
% E field strength: files ending on .E.norm
% normal component of E: files ending on .E.normal
fname_res='.../simu/subject_overlays/lh.test_TMS_1-0001_MagVenture_MC_B70_scalar.central.E.normal'; % UDPATE

% optional: output mesh name for visualizing results using gmsh
fname_msh=''; % UDPATE

% optional: annotation or label file (NOTE: has to be empty for CAT12 results)
%
% freesurfer annotation files are in fs_{subID}/label
% {lh, rh}.aparc.annot (Desikan-Killiany Atlas)
% {lh, rh}.aparc.a2009s.annot (Destrieux Atlas)
% {lh, rh}.aparc.DKTatlas40.annot (DKT Atlas)
% {lh, rh}.BA.annot (some Brodmann areas, unthresholded)
% {lh, rh}.BA.thresh.annot (some Brodmann areas, thresholded)
fname_annot=''; % UDPATE


% read in surfaces
%-------------------------
m=mesh_load_fssurf(fname_surface,bothHemi,true);


% add results as node_data to mesh
%----------------------------------------------
fname_hlp=dir(fname_res);
name_res=fname_hlp.name(find(fname_hlp.name=='.',1)+1:end);

m=mesh_load_fsresults(m, name_res, fname_res ,bothHemi);

if ~isempty(fname_msh)
    mesh_save_gmsh4(m,fname_msh);
end


% get extrema and percentiles for results
%-------------------------------------------------------
if isempty(fname_annot)
    res = mesh_get_surf_extrema_and_percentiles(m,1,[]);
    % m:    input mesh with results data attached
    % 1:    results data is first node data entry
    % []:   use "standard" percentiles
else
    % add annotations or labels as 2nd node_data to mesh
    [m, struct_names]=mesh_load_fsannot(m,fname_annot,bothHemi);
    
    [res, res_annot] = mesh_get_surf_extrema_and_percentiles(m,1,[], 2, struct_names);
    % m:    input mesh with results data attached
    % 1:    results data is first node data entry
    % []:   use "standard" percentiles
    % 2:    annotations are 2nd node data entry
    % struct_names: names of annotations or labels
end


% display results
%-------------------------------------------------------
disp(' ')
disp('RESULTS:')
disp(' ')
for i=1:length(res)
    disp(res(i))
    disp(' ');
end
                                                            
% NOTE: when annotations were used, res_annot contains results
%       similar to those in res for each label listed in struct_names
