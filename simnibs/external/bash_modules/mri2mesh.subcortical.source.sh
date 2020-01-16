# check for input files
for i in {$M2M_DIR/wm.nii.gz,$M2M_DIR/gm.nii.gz,$M2M_DIR/T1fs_conform.nii.gz}; do
 if [ ! -f $i ]; then 
  echo "ERROR: Input file '$i' doesn't exist. Did you run --brain/--brainf?"
  exit 1
 fi
done

OLD_PWD=`pwd`
e2 cd $M2M_DIR/tmp


# ===========================================================
SAY "creating ventricles mask"
# ===========================================================
if $KEEP_MASKS; then
  SAY "using existing ventricles_raw volume!"
  if [ ! -f $MASK_DIR/ventricles_raw.nii.gz ]; then 
    echo "ERROR: Did not find $MASK_DIR/ventricles_raw.nii.gz"
    exit 1
  fi
  e imcp $MASK_DIR/ventricles_raw ventricles_raw

else
  LABELS_VENTRICLES=(4 5 14 24 31 43 44 63)

  # get ventricles mask from FS results
  e mri_convert $FS_DIR/mri/aparc.a2009s+aseg.mgz ./segmentation-mask.nii.gz
  e fslswapdim segmentation-mask.nii.gz x z -y segmentation-mask.nii.gz
  e fslmaths segmentation-mask.nii.gz -mul 0 -bin ventricles_raw.nii.gz
  for i in ${LABELS_VENTRICLES[*]}; do
    e fslmaths segmentation-mask.nii.gz -thr $i -uthr $i -add ventricles_raw.nii.gz -bin ventricles_raw.nii.gz
  done

  # clean mask 
  e fslmaths ventricles_raw.nii.gz -dilM -ero -bin ventricles_raw.nii.gz

  # store in mask_prep
  e imcp ventricles_raw $MASK_DIR/ventricles_raw
fi

# ensure that ventricles are fully inside shrunken WM (volume decoupling)
e fslmaths $M2M_DIR/wm.nii.gz -ero -bin wm_small.nii.gz
e fslmaths ventricles_raw.nii.gz -mul wm_small.nii.gz -bin ventricles_raw.nii.gz

# convert ventricles mask to FS space
e fslswapdim ventricles_raw.nii.gz x -z y ventricles_FS.nii.gz
e mri_convert -odt uchar ./ventricles_FS.nii.gz ./ventricles_FS.nii.gz

# create ventricles surface from volume, smooth surface and convert to stl format
e mri_tessellate ./ventricles_FS.nii.gz 255 ./ventricles_FS.fsmesh
e mris_smooth -n 5 ./ventricles_FS.fsmesh ./ventricles_FS.fsmesh
e mris_convert ./ventricles_FS.fsmesh ./ventricles_FS.stl

# meshfix steps for ventricles surface
e ${BINDIR}/meshfix ventricles_FS.stl -a 2.0 -q -o ventricles.off
e ${BINDIR}/meshfix ventricles.off -a 2.0 -u 5 --vertices $[$NUMBER_OF_VERTICES/8] -q -o ventricles.off
e ${BINDIR}/meshfix ventricles.off -a 2.0 --smooth 1 -q -o ventricles.off
e ${BINDIR}/meshfix ventricles.off -a 2.0 -u 5 -q -o ventricles.off
# next steps included for final cleaning (volume meshing failed in one case otherwise)
e ${BINDIR}/meshfix ventricles.off -a 2.0 -q -o ventricles.off
e ${BINDIR}/meshfix ventricles.off -a 2.0 -u 1 -q -o ventricles.off
e ${BINDIR}/meshfix ventricles.off -a 2.0 -q -o ventricles.off

# decouple ventricles from wm (should not be necessary after volume decoupling; included just for safety)
while ${BINDIR}/meshfix ventricles.off $M2M_DIR/wm.stl --shells 2 --no-clean --intersect; do
 e ${BINDIR}/meshfix ventricles.off $M2M_DIR/wm.stl -a 2.0 --shells 2 --decouple-inin 1 -o ventricles.off
 e ${BINDIR}/meshfix ventricles.off $M2M_DIR/wm.stl -a 2.0 --shells 2 --cut-outer 0 -o ventricles.off
 e ${BINDIR}/meshfix ventricles.off -a 2.0 -u 1 -o ventricles.off
done

e ${BINDIR}/meshfix ventricles.off --stl -q -o $M2M_DIR/ventricles.stl


# ===========================================================
SAY "creating brainstem and cerebellum mask"
# ===========================================================
if $KEEP_MASKS; then
  SAY "using existing cerebellum_raw volume!"
  if [ ! -f $MASK_DIR/cerebellum_raw.nii.gz ]; then 
    echo "ERROR: Did not find $MASK_DIR/cerebellum_raw.nii.gz"
    exit 1
  fi
  e imcp $MASK_DIR/cerebellum_raw cerebellum_raw

else
  LABELS_CEREBELLUM=(7 8 16 46 47)

  # get cerebellum mask from FS results
  e fslmaths segmentation-mask.nii.gz -mul 0 -bin cerebellum_raw.nii.gz
  for i in ${LABELS_CEREBELLUM[*]}; do
    e fslmaths segmentation-mask.nii.gz -thr $i -uthr $i -add cerebellum_raw.nii.gz -bin cerebellum_raw.nii.gz
  done
  e fslmaths cerebellum_raw -bin cerebellum_raw

  # store in mask_prep
  e imcp cerebellum_raw $MASK_DIR/cerebellum_raw
fi

# ensure that mask is fully outside enlarged GM (volume decoupling)
e fslmaths $M2M_DIR/gm.nii.gz -dilM -bin gm_large.nii.gz
e fslmaths cerebellum_raw.nii.gz -sub gm_large.nii.gz -bin cerebellum_raw.nii.gz
e fslmaths cerebellum_raw.nii.gz -ero -dilM -bin cerebellum_raw.nii.gz

# convert cerebellum mask to FS space and format
e fslswapdim cerebellum_raw.nii.gz x -z y cerebellum_FS.nii.gz
e mri_convert -odt uchar ./cerebellum_FS.nii.gz ./cerebellum_FS.nii.gz

# create cerebellum surface from volume, smooth surface and convert to stl
e mri_tessellate ./cerebellum_FS.nii.gz 255 ./cerebellum_FS.fsmesh
e mris_smooth -n 5 ./cerebellum_FS.fsmesh ./cerebellum_FS.fsmesh
e mris_convert ./cerebellum_FS.fsmesh ./cerebellum_FS.stl

# meshfix steps for cerebellum surface
e ${BINDIR}/meshfix cerebellum_FS.stl --no-clean -q -o cerebellum.off
e ${BINDIR}/meshfix cerebellum.off -a 2.0 -q -o cerebellum.off
e ${BINDIR}/meshfix cerebellum.off -a 2.0 -u 5 --vertices $[$NUMBER_OF_VERTICES/7] -q -o cerebellum.off
e ${BINDIR}/meshfix cerebellum.off -a 2.0 -q -o cerebellum.off # these steps included for final cleaning
e ${BINDIR}/meshfix cerebellum.off -a 2.0 -u 1 -q -o cerebellum.off
e ${BINDIR}/meshfix cerebellum.off -a 2.0 -q -o cerebellum.off

# decouple cerebellum from gm (should not be necessary after volume decoupling; included just for safety)
while ${BINDIR}/meshfix cerebellum.off gm.stl --shells 2 --no-clean --intersect; do
 e ${BINDIR}/meshfix cerebellum.off gm.stl -a 2.0 --shells 2 --decouple-outin 0 -o cerebellum.off
 e ${BINDIR}/meshfix cerebellum.off gm.stl -a 2.0 --shells 2 --cut-inner 0 -o cerebellum.off
 e ${BINDIR}/meshfix cerebellum.off -a 2.0 -u 1 -o cerebellum.off
done

e ${BINDIR}/meshfix cerebellum.off --stl -q -o $M2M_DIR/cerebellum.stl


# ===========================================================
SAY "creating volume masks from final ventricles and cerebellum surfaces"
# ===========================================================
e ${BINDIR}/meshfix $M2M_DIR/ventricles.stl --no-clean --fsmesh -o ventricles.fsmesh
e mris_transform --dst ../ref_FS.nii.gz --src ../ref_FS.nii.gz ./ventricles.fsmesh ../unity.xfm ./ventricles.fsmesh
e mris_fill -c ./ventricles.fsmesh $M2M_DIR/ventricles.nii.gz
e fslswapdim $M2M_DIR/ventricles.nii.gz x z -y $M2M_DIR/ventricles.nii.gz

e ${BINDIR}/meshfix $M2M_DIR/cerebellum.stl --no-clean --fsmesh -o cerebellum.fsmesh
e mris_transform --dst ../ref_FS.nii.gz --src ../ref_FS.nii.gz ./cerebellum.fsmesh ../unity.xfm ./cerebellum.fsmesh
e mris_fill -c ./cerebellum.fsmesh $M2M_DIR/cerebellum.nii.gz
e fslswapdim $M2M_DIR/cerebellum.nii.gz x z -y $M2M_DIR/cerebellum.nii.gz

e2 cd $OLD_PWD

