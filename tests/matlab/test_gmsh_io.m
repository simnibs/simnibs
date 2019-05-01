% test of:
%  mesh_load_gmsh4
%  mesh_save_gmsh4
%
% AT 09-Apr-2018

clear all
close all

testmsh='sphere.msh';
testbinmsh='sphere_bin.msh';
testmsh4='sphere_gmsh4.msh';
testbinmsh4='sphere_gmsh4_bin.msh';
testdir=fileparts(mfilename('fullpath'));

testres=[];

% test the skipping of elements that are neither
% triangles nor tetrahedra
disp('---------------------------------------------------------')
disp('Test: Reading ascii and binary mesh with points and lines')
disp('Note: This gives warnings that some elements are skipped')
disp('---------------------------------------------------------')

m=mesh_load_gmsh4([testdir filesep testmsh]);
m2=mesh_load_gmsh4([testdir filesep testbinmsh]);


if anydifferences(m,m2)
    disp('------- test NOT PASSED ------');
    testres=[testres false];
else
    disp('------- test OK ------');
    testres=[testres true];
end

m3=mesh_load_gmsh4([testdir filesep testmsh4]);

% if anydifferences(m,m3)
%     disp('------- test NOT PASSED ------');
%     testres=[testres false];
% else
%     disp('------- test OK ------');
%     testres=[testres true];
% end

m3=mesh_load_gmsh4([testdir filesep testbinmsh4]);
% if anydifferences(m,m3)
%     disp('------- test NOT PASSED ------');
%     testres=[testres false];
% else
%     disp('------- test OK ------');
%     testres=[testres true];
% end

% create node data and element data, write (ascii and binary), read and compare
disp('---------------------------------------------------------')
disp('Test: Reading and writing of node and element data')
disp('---------------------------------------------------------')

m.node_data{1}.name='test';
m.node_data{1}.data=rand(size(m.nodes));
m.node_data{2}.name='test2';
m.node_data{2}.data=rand(size(m.nodes,1),1);

m.element_data{1}.name='test3';
m.element_data{1}.tridata=rand(size(m.triangles));
m.element_data{1}.tetdata=rand(size(m.tetrahedra,1),3);

m.element_data{2}.name='test4';
m.element_data{2}.tridata=rand(size(m.triangles,1),1);
m.element_data{2}.tetdata=rand(size(m.tetrahedra,1),1);

disp('testing binary')
mesh_save_gmsh4(m,[testdir filesep 'test123']);
m2=mesh_load_gmsh4([testdir filesep 'test123.msh']);

if anydifferences(m,m2)
    disp('------- test NOT PASSED ------');
    testres=[testres false];
else
    disp('------- test OK ------');
    testres=[testres true];
end

disp('testing ascii')
mesh_save_gmsh4(m,[testdir filesep 'test123'],'ascii');
m2=mesh_load_gmsh4([testdir filesep 'test123.msh']);

if anydifferences(m,m2)
    disp('------- test NOT PASSED ------');
    testres=[testres false];
else
    disp('------- test OK ------');
    testres=[testres true];
end;


% write only tet element data, write (ascii and binary), read and compare
disp('---------------------------------------------------------')
disp('Test: Reading and writing of node and element data (tetdata only)')
disp('---------------------------------------------------------')

m.element_data{1}.name='test3';
m.element_data{1}.tridata=[];
m.element_data{1}.tetdata=rand(size(m.tetrahedra,1),3);

m.element_data{2}.name='test4';
m.element_data{2}.tridata=[];
m.element_data{2}.tetdata=rand(size(m.tetrahedra,1),1);

disp('testing binary')
mesh_save_gmsh4(m,[testdir filesep 'test123']);
m2=mesh_load_gmsh4([testdir filesep 'test123.msh']);

if anydifferences(m,m2)
    disp('------- test NOT PASSED ------');
    testres=[testres false];
else
    disp('------- test OK ------');
    testres=[testres true];
end;

disp('testing ascii')
mesh_save_gmsh4(m,[testdir filesep 'test123'],'ascii');
m2=mesh_load_gmsh4([testdir filesep 'test123.msh']);

if anydifferences(m,m2)
    disp('------- test NOT PASSED ------');
    testres=[testres false];
else
    disp('------- test OK ------');
    testres=[testres true];
end;


% write only tri element data, write (ascii and binary), read and compare
disp('---------------------------------------------------------')
disp('Test: Reading and writing of node and element data (tridata only)')
disp('---------------------------------------------------------')

m.element_data{1}.name='test3';
m.element_data{1}.tridata=rand(size(m.triangles,1),3);
m.element_data{1}.tetdata=[];

m.element_data{2}.name='test4';
m.element_data{2}.tridata=rand(size(m.triangles,1),1);
m.element_data{2}.tetdata=[];

disp('testing binary')
mesh_save_gmsh4(m,[testdir filesep 'test123']);
m2=mesh_load_gmsh4([testdir filesep 'test123.msh']);

if anydifferences(m,m2)
    disp('------- test NOT PASSED ------');
    testres=[testres false];
else
    disp('------- test OK ------');
    testres=[testres true];
end;

disp('testing ascii')
mesh_save_gmsh4(m,[testdir filesep 'test123'],'ascii');
m2=mesh_load_gmsh4([testdir filesep 'test123.msh']);

if anydifferences(m,m2)
    disp('------- test NOT PASSED ------');
    testres=[testres false];
else
    disp('------- test OK ------');
    testres=[testres true];
end;



% write a mesh with tetrahedra only
disp('---------------------------------------------------------')
disp('Test: Reading and writing of a mesh with tetrahedra only')
disp('---------------------------------------------------------')

m.element_data{1}.tetdata=rand(size(m.tetrahedra,1),3); % add some tet element data again
m.element_data{2}.tetdata=rand(size(m.tetrahedra,1),1);

m2=mesh_extract_regions(m,'tet');

disp('testing binary')
mesh_save_gmsh4(m2,[testdir filesep 'test123']);
m3=mesh_load_gmsh4([testdir filesep 'test123.msh']);

if anydifferences(m2,m3)
    disp('------- test NOT PASSED ------');
    testres=[testres false];
else
    disp('------- test OK ------');
    testres=[testres true];
end;

disp('testing ascii')
mesh_save_gmsh4(m2,[testdir filesep 'test123'],'ascii');
m3=mesh_load_gmsh4([testdir filesep 'test123.msh']);

if anydifferences(m2,m3)
    disp('------- test NOT PASSED ------');
    testres=[testres false];
else
    disp('------- test OK ------');
    testres=[testres true];
end;



% write a mesh with triangles only
disp('---------------------------------------------------------')
disp('Test: Reading and writing of a mesh with tetrahedra only')
disp('---------------------------------------------------------')

m2=mesh_extract_regions(m,'tri');

disp('testing binary')
mesh_save_gmsh4(m2,[testdir filesep 'test123']);
m3=mesh_load_gmsh4([testdir filesep 'test123.msh']);

if anydifferences(m2,m3)
    disp('------- test NOT PASSED ------');
    testres=[testres false];
else
    disp('------- test OK ------');
    testres=[testres true];
end;

disp('testing ascii')
mesh_save_gmsh4(m2,[testdir filesep 'test123'],'ascii');
m3=mesh_load_gmsh4([testdir filesep 'test123.msh']);

if anydifferences(m2,m3)
    disp('------- test NOT PASSED ------');
    testres=[testres false];
else
    disp('------- test OK ------');
    testres=[testres true];
end;

delete([testdir filesep 'test123.msh'])
cd(testdir)

disp(' ');
if all(testres)

    disp('------- ALL TESTS OK ------');
else

    error('------- at lest one test not passed ------');
end;

