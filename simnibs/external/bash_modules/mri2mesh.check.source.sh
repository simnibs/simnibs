# check for the input files
for i in {{wm,gm,cerebellum,csf,skull,skin,ventricles}.stl,T1fs_conform.nii.gz}; do
if [ ! -f $M2M_DIR/$i ]; then echo "ERROR: File $M2M_DIR/$i doesn't exist yet. Did you run mri2mesh --all?"; exit 1; fi;
done


if [ $QUICKCHECK = 1 ]; then
  OLD_PWD=`pwd`
  e2 cd $M2M_DIR
  e freeview -v skin.nii.gz:opacity=0.4:colormap=jet:colorscale=0,1 skull.nii.gz:opacity=0.6:colormap=jet:colorscale=0,1.3   csf.nii.gz:opacity=0.6:colormap=jet:colorscale=0,200  wm.nii.gz:opacity=0.7:colormap=jet:colorscale=0,1.3  gm.nii.gz:opacity=0.5:colormap=jet:colorscale=0,2  cerebellum.nii.gz:opacity=0.6:colormap=jet:colorscale=0,2  ventricles.nii.gz:opacity=0.6:colormap=jet:colorscale=0,200 T1fs_conform.nii.gz:opacity=0.6&
   # note: replaced fslview by freeview, as fslview is deprecated.
   # colorscale max corresponds to: 1 red, 1.3 yellow, 200 blue, 2 green

  if [ -d $M2M_DIR/toMNI ]; then 
    e freeview -v $STD_MNI $M2M_DIR/toMNI/T1fs_nu_12DOF_MNI.nii.gz:visible=0 $M2M_DIR/toMNI/T1fs_nu_nonlin_MNI.nii.gz &
  fi
  e2 cd $OLD_PWD
  return;
fi;


if [ ! -f $M2M_DIR/tmp/check_tmp1.fsmesh -o $check_cache -eq 0 ]; then
 SAY "Merging stl files and converting to freesurfer mesh format"
 e ${BINDIR}/meshfix $M2M_DIR/wm.stl $M2M_DIR/gm.stl --shells 2 --no-clean -o $M2M_DIR/tmp/check_tmp1.off
 e ${BINDIR}/meshfix $M2M_DIR/tmp/check_tmp1.off $M2M_DIR/cerebellum.stl --shells 3 --no-clean -o $M2M_DIR/tmp/check_tmp1.off
 e ${BINDIR}/meshfix $M2M_DIR/tmp/check_tmp1.off $M2M_DIR/ventricles.stl --shells 5 --no-clean -o $M2M_DIR/tmp/check_tmp1.off
 e ${BINDIR}/meshfix $M2M_DIR/tmp/check_tmp1.off $M2M_DIR/csf.stl --shells 6 --no-clean -o $M2M_DIR/tmp/check_tmp1.off
 e ${BINDIR}/meshfix $M2M_DIR/tmp/check_tmp1.off $M2M_DIR/skull.stl --shells 7 --no-clean -o $M2M_DIR/tmp/check_tmp1.off
 e ${BINDIR}/meshfix $M2M_DIR/tmp/check_tmp1.off $M2M_DIR/skin.stl --shells 8 --no-clean -o $M2M_DIR/tmp/check_tmp1.off

 e ${BINDIR}/meshfix $M2M_DIR/tmp/check_tmp1.off  --no-clean --shells 8 --fsmesh -o $M2M_DIR/tmp/check_tmp1.fsmesh
 e mris_transform --dst $M2M_DIR/ref_FS.nii.gz --src $M2M_DIR/ref_FS.nii.gz $M2M_DIR/tmp/check_tmp1.fsmesh $M2M_DIR/unity.xfm $M2M_DIR/tmp/check_tmp1.fsmesh
fi

VOLUMEFILES=$M2M_DIR/'T1fs_nu_conform.nii.gz:grayscale=50,230'
if [ -f $M2M_DIR/T1_conform.nii.gz ]; then
  VOLUMEFILES=`echo $M2M_DIR/T1_conform.nii.gz $VOLUMEFILES `
fi
if [ -f $M2M_DIR/T2_conform.nii.gz ]; then
  VOLUMEFILES=`echo $M2M_DIR/T2_conform.nii.gz $VOLUMEFILES `
fi
if [ -f $M2M_DIR/T2fs_conform.nii.gz ]; then
  VOLUMEFILES=`echo $M2M_DIR/T2fs_conform.nii.gz $VOLUMEFILES `
fi

e freeview -v $VOLUMEFILES -f $FS_DIR/surf/lh.pial:edgethickness=1 $FS_DIR/surf/lh.white:edgethickness=1 $FS_DIR/surf/rh.pial:edgethickness=1 $FS_DIR/surf/rh.white:edgethickness=1 $M2M_DIR/tmp/check_tmp1.fsmesh:edgecolor=red & 

if [ -d $M2M_DIR/toMNI ]; then 
  e freeview -v $STD_MNI $M2M_DIR/toMNI/T1fs_nu_12DOF_MNI.nii.gz:visible=0 $M2M_DIR/toMNI/T1fs_nu_nonlin_MNI.nii.gz &
fi










