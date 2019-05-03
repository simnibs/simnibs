%
% Example script that reads a SimNIBS mesh with simulation results and
% extracts the gray matter tetrahedra.
%
% To characterize the reached field strength, it determines the
% peak magnitude (which might be affected by outliers) and some percentiles
% of normE.
% 
% To characterize focality, it determines the GM volue in which the field
% exceeds a certain fraction of the 99.9 percentile.
%
% Finally, it plots a histogram of the field values and the field on the
% GM surface
%
% A. Thielscher, 11-Apr-2018; updated 13-Sep-2018


% load mesh
[fname,pname] = uigetfile('*.msh','Select a mesh with simulation results');
if isequal(fname,0) || isequal(pname,0); return; end
m_head=mesh_load_gmsh4([pname fname]);

% extract GM tetrahedra (GM has region number 2) and their data
% and find data field containing the magnitude of E
m=mesh_extract_regions(m_head,'tet',2);

idxNormE=0;
for i=1:length(m.element_data)
    if strcmpi(m.element_data{i}.name,'normE')
        idxNormE=i;
        break;
    end
end


% measures of field intensity
%--------------------
% get maximum field strength in GM
normE=m.element_data{idxNormE}.tetdata;
disp(' ');
disp('---------------------------------- ');
disp(['The peak field strength in GM is ' num2str(max(normE)) ' V/m']);

% get some percentiles of the field strength:
% get tetrahedra volumes and sort the volumes using the 
% field strength values to get the percentiles
percentiles=[50 75 90 95 99 99.5 99.9];

tetvols=mesh_get_tetrahedron_sizes(m);

[normE,idx] = sort(normE);
tetvols=tetvols(idx);
tetvols=cumsum(tetvols);
tetvolsNorm=tetvols/tetvols(end); % normalize to one

for i=1:length(percentiles)
    idx=find(tetvolsNorm>percentiles(i)/100,1,'first');
    normE_prc(i)=normE(idx);
end

disp(' ');
disp('percentiles of normE in GM in V/m:')
disp(num2str(percentiles));
disp(num2str(normE_prc));


% measure of focality
%--------------------
% GM volume in which the field exceeds 50%, 75% and 90% of the 99.9
% percentile
cutoffs=[50 75 90];
for i=1:length(cutoffs)
    idx=find(normE>cutoffs(i)/100*normE_prc(end),1,'first');
    GMvol_cutoff(i)=tetvols(end)-tetvols(idx);
end

disp(' ');
disp('GM volume (in mm³), in which normE exceeds the');
disp(['stated fractions of the ' num2str(percentiles(end)) ' percentile']);
disp(num2str(cutoffs));
disp(num2str(GMvol_cutoff,'%.0f'));
disp('---------------------------------------- ');
disp(' ');

% plot histogram of field in GM volume and show field on GM surface
%--------------------
figure;
h = subplot(1,2,1);
mesh_get_histogram(m,'haxis',h)
h = subplot(1,2,2);
mesh_show_surface(m_head,'haxis',h)
