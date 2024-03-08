function test_gifti

N = 31;
[x, y] = meshgrid (1:N);
tri = delaunay (x(:), y(:));
z = peaks (N);
p = struct ('faces', tri, 'vertices', [x(:) y(:) z(:)]);

g = gifti(p);

basename = tempname;
file = [basename '.gii'];
save(g,file,'ASCII');
g = gifti(file);
delete(file);
save(g,file,'Base64Binary');
g = gifti(file);
delete(file);
save(g,file,'GZipBase64Binary');
g = gifti(file);
delete(file);
save(g,file,'ExternalFileBinary');
g = gifti(file);
delete(file);
delete([basename '.dat']);
