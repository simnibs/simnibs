function notpass = anydifferences(m,m2)
% helper function  for test_gmsh_io.m

notpass = false;

if max(abs(m.nodes(:)-m2.nodes(:)))>1e-6
    warning(['max. difference in node positions: ' num2str(max(abs(m.nodes(:)-m2.nodes(:))))])
    notpass=true;
end;

if any(m.triangles(:) ~= m2.triangles(:));
    warning('triangle indices changed!');
    notpass=true; return;
end;

if any(m.triangle_regions(:) ~= m2.triangle_regions(:));
    warning('triangle regions changed!');
    notpass=true; return;
end;

if any(m.tetrahedra(:) ~= m2.tetrahedra(:));
    warning('tetrahedra indices changed!');
    notpass=true; return;
end;

if any(m.tetrahedron_regions(:) ~= m2.tetrahedron_regions(:));
    warning('tetrahedra indices changed!');
    notpass=true; return;
end;

for i=1:length(m.node_data)
    if ~strcmp(m.node_data{i}.name,m2.node_data{i}.name)
        warning('name of node data fields does not match');
        notpass=true; return; 
    end;
    if max(abs(m.node_data{i}.data(:)-m2.node_data{i}.data(:)))>eps(single(1))
        warning(['max. difference in node data ' num2str(i) ' ' ...
                 num2str(max(abs(m.node_data{i}.data(:)-m2.node_data{i}.data(:))))])
        notpass=true;
    end;
end;

for i=1:length(m.element_data)
    if ~strcmp(m.element_data{i}.name,m2.element_data{i}.name)
        warning('name of element data fields does not match');
        notpass=true; return;
    end;
    
    if max(abs(m.element_data{i}.tridata(:)-m2.element_data{i}.tridata(:)))>eps(single(1))
        warning(['max. difference in element tridata ' num2str(i) ' ' ...
                 num2str(max(abs(m.element_data{i}.tridata(:)-m2.element_data{i}.tridata(:))))])
        notpass=true;
    end;
    
    if max(abs(m.element_data{i}.tetdata(:)-m2.element_data{i}.tetdata(:)))>eps(single(1))
        warning(['max. difference in element tetdata ' num2str(i) ' ' ...
                 num2str(max(abs(m.element_data{i}.tetdata(:)-m2.element_data{i}.tetdata(:))))])
        notpass=true; 
    end;   
end;

