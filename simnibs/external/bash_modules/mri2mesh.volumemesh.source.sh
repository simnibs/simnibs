# check for first input file
if [ ! -d $M2M_DIR -o ! -f $M2M_DIR/wm.stl ]; then echo "ERROR: Input file doesn't exist. Did you call this script from mri2mesh?"; exit 1; fi;

e2 cat > $SUBJECT.geo <<EOF
//to mesh the volumes call ${BINDIR}/gmsh with
//gmsh -3 -o $SUBJECT.msh $SUBJECT.geo
Mesh.Algorithm3D=4; //1=delaunay (tetgen) and 4=frontal (netgen)
Mesh.Optimize=1;
Mesh.OptimizeNetgen=1;

Merge "$M2M_DIR/wm.stl"; // first surface
Merge "$M2M_DIR/gm.stl";
Merge "$M2M_DIR/csf.stl";
Merge "$M2M_DIR/skull.stl";
Merge "$M2M_DIR/skin.stl";
Merge "$M2M_DIR/cerebellum.stl";
Merge "$M2M_DIR/ventricles.stl";

Surface Loop(1) = {1}; // surface number on rhs; 1: wm.stl
Surface Loop(2) = {2}; // gm.stl
Surface Loop(3) = {3}; // csf.stl
Surface Loop(4) = {4}; // skull.stl
Surface Loop(5) = {5}; // skin.stl
Surface Loop(6) = {6}; // cerebellum.stl
Surface Loop(7) = {7}; // ventricles.stl

Volume(1) = {1, 7};    // surface loop numbers on rhs; e.g., volume between surface loops 1 (wm) and 7 (ventricles) is sm
Volume(2) = {1, 2};    // gm volume
Volume(3) = {2, 6, 3}; // csf volume: outside gm, cerebellum; inside csf
Volume(4) = {3, 4};    // skull volume
Volume(5) = {4, 5};    // skin volume
Volume(6) = {6};       // cerebellum volume
Volume(7) = {7};       // ventricle volume

Physical Surface(1001) = {1,6}; // merge wm and cerebellum; LHS: target surface region number, RHS: surface number (i.e. from merge ...)
Physical Surface(1002) = {2}; 
Physical Surface(1003) = {3,7}; // merge csf and ventricles
Physical Surface(1004) = {4};
Physical Surface(1005) = {5};

Physical Volume(1) = {1,6}; // merge wm and cerebellum; LHS: target volume region number, RHS: volume number
Physical Volume(2) = {2}; 
Physical Volume(3) = {3,7}; // merge csf and ventricles
Physical Volume(4) = {4};
Physical Volume(5) = {5};

EOF

e ${BINDIR}/gmsh -3 -bin -o $SUBJECT.msh $SUBJECT.geo

# fix labeling, thin tetrahedra, normal directions
e merge_labels_WM $SUBJECT'.msh' $SUBJECT'.msh'

# create final wm and gm masks from mesh, using the same resampling as for the later FEM results
e msh2nii --create_masks $SUBJECT'.msh' $M2M_DIR/T1fs_nu_conform.nii.gz $M2M_DIR/tmp/MasksfromMesh
e immv $M2M_DIR/tmp/MasksfromMesh_mask_1 $M2M_DIR/wm_fromMesh
e immv $M2M_DIR/tmp/MasksfromMesh_mask_2 $M2M_DIR/gm_fromMesh

# Final check of the mesh
e check_meshed_volumes -m $M2M_DIR --remove
