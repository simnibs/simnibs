function [S,SH] = cat_surf_GI3D(S,D,R,opt)
% =========================================================================
% Next version of local gyrification index calculation. The laplace-method
% is used to map the area of the face of a given surface S to the hull
% given by D and R. Because lapace mapping the vertices with equal
% triangulation produce errors (at watershed of laplace field in gyral
% regions where on e vertex move to the other side and the face get
% corrupt.) we only map vertices to their hull position and retriangulate
% the surface with delaunay. Next, voroi area of this vertices is estimated
% to get the outer surface area (OSA) and the inner surface area (ISA) to
% calculate the GI.
% _________________________________________________________________________
% INPUT:
%   S.faces    = Surface faces
%   S.vertices = Surface points
%   D          = Tissues image for laplace-filter
%   R          = filter region 
%   opt        = options...
%     .          = 0.02 (lapalce error criteria)
%     .streamopt = [,] (stepsize and maxiumum length of a streamline)
%
% OUTPUT:
%   S.faces       = Surface faces
%   S.vertices    = Surface points
%   S.area        = Surface area for each vertex
%   SH.faces      = Hull surface faces
%   SH.vertices   = Hull surface faces
%   SH.area       = Hull surface area for each vertex
%
%   S.GI          = Gyrifcation Index for points (S.area / S.Harea)
%   S.SL          = streamlength between S.vertices and S.Hvertices
% =========================================================================

  % pinterpol  = (2^opt.interpol)-(2^opt.interpol-1)*(opt.interpol>0);
  % if ~isfield(opt,'gridres'), opt.gridres = [1,1,1]; end
  
  % ATTENSION matlab-stream function only work with double! so don't
  % convert vertices!
  S.faces = single(S.faces);
   
  D = single(D);
  R = single(R);
  
  % for correct smoother border (else WM blows to strong
  D = D - 0.5*(single(D==1) & imdilate(D==0,ones(3,3)));
  
  % streamline calculation 
  % _______________________________________________________________________
  if ~exist('opt','var') opt = struct(); end
  def.side          = 'streams0';
  def.interpV       = 1; 
  def.streamopt(1)  = 0.01; % 0.05 
  def.streamopt(2)  = 1000 ./ def.streamopt(1);
  opt = cat_io_checkinopt(opt,def); 
  
  SL = cat_surf_epivolsurf(D,R,opt,S);
  SL.OP = SL.L(:,:,7);
  S.Hvertices  = SL.OP;
  S.SL         = SL.SL; 
  S.streams    = SL.streams;
  
  
  SH.vertices = SL.OP;
  SH.faces    = S.faces;
  SH.facevertexcdata = SL.SL;
  
  % GI: gyrification index 
  % _______________________________________________________________________
  % Calculate the face areas of the surface and it's hull (see surfacearea)
  % and map them to the vertices (see verticemapping).
  % The GI is the ration between the inner and outer area.
  S = surfacearea(S,'');  S.farea  = S.area;  S = verticemapping(S,'area'); 
  S = surfacearea(S,'H'); S.Hfarea = S.Harea; S = verticemapping(S,'Harea');
  S.GI  = S.area ./ S.Harea; 
  S.GIl = log(S.GI);
  %fprintf(1,'|GI:%0.0f',toc); tic

end


function S = surfacearea(S,fname)
% calculate surface area of the faces and also map it to the vertices
% IN:   S.*vertices, S.faces
% OUT:  S.*verticesarea, S.*facearea

  fa = [fname 'area']; fd = [fname 'dist'];
  ndim = size(S.vertices,2);
  
  % facearea (Horonsche Form)
  S = surfacedistance(S,fname);
  if ndim==2
    S.(fa) = S.(fd);
  elseif ndim==3
    facesp = sum(S.(fd),2) / 2;  % s = (a + b + c) / 2;
    S.(fa) = (facesp .* (facesp - S.(fd)(:,1)) .* (facesp - S.(fd)(:,2)) .* (facesp - S.(fd)(:,3))).^0.5; % area=sqrt(s*(s-a)*(s-b)*(s-c));
  end
  % numberical (to small point diffences) and mapping broblems (crossing of streamlines)
  % -> correction because this is theoretical not possible (laplace field theorie)
  S.(fa)(S.(fa)==0) = eps; % to small values
  S.(fa) = abs(S.(fa));    % streamline-crossing
end
function S = verticemapping(S,fname)
% mapping of facedata to vertices

  data  = zeros(size(S.vertices,1),1);
  [v,f] = sort(S.faces(:)); [f,fj] = ind2sub(size(S.faces),f);  %#ok<NASGU>
  far = S.(fname)(f);
  for i=1:numel(S.faces), data(v(i)) = data(v(i)) + far(i); end

%   for i=1:numel(S.faces), S.(va)(v(i)) = S.(va)(v(i)) + S.(fa)(f(i)); end
  
%   if ndim==2
%     for i=1:size(S.vertices,1)
%       p=S.faces(i,1); S.(va)(p) = S.(va)(p) + S.(fa)(i);
%       p=S.faces(i,2); S.(va)(p) = S.(va)(p) + S.(fa)(i);
%     end
%   elseif ndim==3
%     for i=1:size(S.faces,1)
%       p=S.faces(i,1); S.(va)(p) = S.(va)(p) + S.(fa)(i);
%       p=S.faces(i,2); S.(va)(p) = S.(va)(p) + S.(fa)(i);
%       p=S.faces(i,3); S.(va)(p) = S.(va)(p) + S.(fa)(i);
%     end
%   end
%S.(va) = arrayfun(@(p,i) S.(va)(p) + S.(fa)(i),S.faces,(1:size(S.faces,1))'*ones(1,ndim));

  data = data / size(S.vertices,2); % Schwerpunkt... besser Voronoi, aber wie bei ner Oberfläche im Raum???
  S = rmfield(S,fname);
  S.(fname) = data;
end

%{
function S = verticeneighbor(S)
  if ~isfield(S,'nb')
    S.nb = [S.vertices(:,1),S.vertices(:,2);
            S.vertices(:,2),S.vertices(:,3);
            S.vertices(:,3),S.vertices(:,1)];
  end
end
%}

function S = surfacedistance(S,fname)
% 2D: d(AB)
% 3D: [c,a,b] = d(AB),d(BC),d(CA)
  v = [fname 'vertices']; fd = [fname 'dist'];
  ndim = size(S.vertices,2);
  
  % facearea (Horonsche Form)
  if ndim==2
    S.(fd) = sum( (S.(v)(S.faces(:,1),:) - S.(v)(S.faces(:,2),:)).^2 , 2) .^ 0.5;
  elseif ndim==3
    S.(fd) = [sum( (S.(v)(S.faces(:,1),:) - S.(v)(S.faces(:,2),:)).^2 , 2) .^ 0.5, ...
              sum( (S.(v)(S.faces(:,2),:) - S.(v)(S.faces(:,3),:)).^2 , 2) .^ 0.5, ...
              sum( (S.(v)(S.faces(:,3),:) - S.(v)(S.faces(:,1),:)).^2 , 2) .^ 0.5]; % [c,a,b] = d(AB),d(BC),d(CA)
  end
end
