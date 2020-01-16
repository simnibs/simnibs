# check for input files
for i in {$M2M_DIR/T1fs_conform.nii.gz,$M2M_DIR/T1fs_roi.nii.gz,$M2M_DIR/cerebellum.nii.gz}; do
 if [ ! -f $i ]; then 
  echo "ERROR: Input file '$i' doesn't exist. Did you run --brain/--brainf and --subcort?"
  exit 1
 fi
done

OLD_PWD=`pwd`
e2 cd $M2M_DIR

# ==================================
SAY "Coregistering images to T1fs_conform"
# ==================================

# create intensity normalized version of T1fs_roi
e fslswapdim T1fs_roi x -z y tmp/T1fs_roi_FS
e mri_nu_correct.mni --n 1 --proto-iters 1000 --distance 50 --no-rescale --i tmp/T1fs_roi_FS.nii.gz --o tmp/T1fs_nu_roi.nii.gz
e fslswapdim tmp/T1fs_nu_roi x z -y tmp/T1fs_nu_roi

# transform all images into conformed space, cut neck
e fslmaths T1fs_roi -bin -s 1 tmp/T1fs_roi-mask
for i in {T2,T2fs,T1}; do
  if [ -f $i.nii.gz ]; then
    e flirt -in $i -ref T1fs_roi -refweight tmp/T1fs_roi-mask -interp sinc -cost mutualinfo -searchcost mutualinfo -dof 6 -omat tmp/$i'2Conform.mat' -out $i'_conform'
    e fslmaths $i'_conform' -mul tmp/T1fs_roi-mask $i'_conform'
  fi
done

# get linear transformations conformed space <-> MNI space
e flirt -dof 12 -in $STD_MNI -ref T1fs_nu_conform -omat tmp/MNI2conform_forbet.mat
e convert_xfm -omat tmp/Conform2MNI.mat -inverse tmp/MNI2conform_forbet.mat


e2 cd $M2M_DIR/tmp

# ==================================
SAY "Get CSF mask from T1fs and T2fs"
# ==================================
if $KEEP_MASKS; then
  SAY "using existing csf_raw volume!"
  if [ ! -f $MASK_DIR/csf_raw.nii.gz ]; then 
    echo "ERROR: Did not find $MASK_DIR/csf_raw.nii.gz"
    exit 1
  fi
  e imcp $MASK_DIR/csf_raw csf_raw

else
  # get initial surface (.vtk) from intensity normalized T1fs_roi image
  e bet T1fs_nu_roi T1fs_conform_betted -e -m -g -0.2 # trying to include sinuses by setting -g -0.2 

  # create surfaces (produces masks, and off surfaces)
  if [ -f $M2M_DIR/T2fs.nii.gz ]; then
    e betsurf -m -o $M2M_DIR/T1fs_roi $M2M_DIR/T2fs_conform T1fs_conform_betted_mesh.vtk Conform2MNI.mat csf_raw
  elif [ -f $M2M_DIR/T2.nii.gz ]; then
    e betsurf -m -o $M2M_DIR/T1fs_roi $M2M_DIR/T2_conform T1fs_conform_betted_mesh.vtk Conform2MNI.mat csf_raw
  else
    e betsurf --t1only -m -o $M2M_DIR/T1fs_roi T1fs_conform_betted_mesh.vtk Conform2MNI.mat csf_raw
  fi
  e mv csf_raw_inskull_mask.nii.gz csf_raw.nii.gz

  # store in mask_prep
  e imcp csf_raw $MASK_DIR/csf_raw
fi

# create enlarged cerebellum volume and add it to enlarged gm volume
e fslmaths $M2M_DIR/gm.nii.gz -dilM -bin gm_large.nii.gz
e fslmaths $M2M_DIR/cerebellum.nii.gz -dilM -add gm_large.nii.gz -bin brain_large.nii.gz

# ensure that CSF mask is fully outside brain_large (volume decoupling)
e fslmaths csf_raw.nii.gz -add brain_large.nii.gz -bin csf_raw.nii.gz

# convert CSF mask to FS space
e fslswapdim csf_raw.nii.gz x -z y csf_FS.nii.gz
e mri_convert -odt uchar ./csf_FS.nii.gz ./csf_FS.nii.gz

# create CSF surface from volume
e mri_tessellate ./csf_FS.nii.gz 255 ./csf_FS.fsmesh
# smooth CSF surface
e mris_smooth -n 5 ./csf_FS.fsmesh ./csf_FS.fsmesh
# convert CSF surface to stl format
e mris_convert ./csf_FS.fsmesh ./csf_FS.stl

# meshfix steps for CSF surface
e ${BINDIR}/meshfix csf_FS.stl -a 2.0 -u 5 --vertices $NUMBER_OF_VERTICES -q -o csf.off
e ${BINDIR}/meshfix csf.off -a 2.0 -q -o csf.off # this steps included for final cleaning
e ${BINDIR}/meshfix csf.off -a 2.0 -u 1 -q -o csf.off
e ${BINDIR}/meshfix csf.off -a 2.0 -q -o csf.off

# decouple CSF surface from gm and cerebellum (should not be necessary after volume decoupling; included just for safety)
while ${BINDIR}/meshfix csf.off $M2M_DIR/gm.stl --shells 2 --no-clean --intersect; do
 e ${BINDIR}/meshfix csf.off $M2M_DIR/gm.stl -a 2.0 --shells 2 --decouple-outout 0 -o csf.off
 e ${BINDIR}/meshfix csf.off $M2M_DIR/gm.stl -a 2.0 --shells 2 --cut-inner 0 -o csf.off
 e ${BINDIR}/meshfix csf.off -a 2.0 -u 1 -o csf.off
done
while ${BINDIR}/meshfix csf.off $M2M_DIR/cerebellum.stl --shells 2 --no-clean --intersect; do
 e ${BINDIR}/meshfix csf.off $M2M_DIR/cerebellum.stl -a 2.0 --shells 2 --decouple-outout 0 -o csf.off
 e ${BINDIR}/meshfix csf.off $M2M_DIR/cerebellum.stl -a 2.0 --shells 2 --cut-inner 0 -o csf.off
 e ${BINDIR}/meshfix csf.off -a 2.0 -u 1 -o csf.off
done

# final csf surface in stl format
e ${BINDIR}/meshfix csf.off --stl -q -o $M2M_DIR/csf.stl

# create volume mask from csf.stl
e ${BINDIR}/meshfix $M2M_DIR/csf.stl --no-clean --fsmesh -o csf.fsmesh
e mris_transform --dst ../ref_FS.nii.gz --src ../ref_FS.nii.gz ./csf.fsmesh ../unity.xfm ./csf.fsmesh
e mris_fill -r 0.5 -c ./csf.fsmesh $M2M_DIR/csf.nii.gz
e fslswapdim $M2M_DIR/csf.nii.gz x z -y $M2M_DIR/csf.nii.gz


# ==================================
SAY "Get skull mask from T1 and T2"
# ==================================
if $KEEP_MASKS; then
  SAY "using existing skull_raw volume!"
  if [ ! -f $MASK_DIR/skull_raw.nii.gz ]; then 
    echo "ERROR: Did not find $MASK_DIR/skull_raw.nii.gz"
    exit 1
  fi
  e imcp $MASK_DIR/skull_raw skull_raw

else
  # create surfaces (produces masks, and off surfaces)
  if [ -f $M2M_DIR/T1.nii.gz ]; then
    if [ -f $M2M_DIR/T2.nii.gz ]; then
      e betsurf -m -o $M2M_DIR/T1_conform $M2M_DIR/T2_conform T1fs_conform_betted_mesh.vtk Conform2MNI.mat skull_raw
    else
      e betsurf --t1only -m -o $M2M_DIR/T1_conform T1fs_conform_betted_mesh.vtk Conform2MNI.mat skull_raw
    fi
    e mv skull_raw_outskull_mask.nii.gz skull_raw.nii.gz
  else
    e mv csf_raw_outskull_mask.nii.gz skull_raw.nii.gz
  fi

  # store in mask_prep
  e imcp skull_raw $MASK_DIR/skull_raw
fi

# create enlarged csf volume
e fslmaths $M2M_DIR/csf.nii.gz -dilM -bin csf_large.nii.gz

# ensure that skull mask is fully outside csf_large (volume decoupling)
e fslmaths skull_raw.nii.gz -add csf_large.nii.gz -bin skull_raw.nii.gz

# convert skull mask to FS space
e fslswapdim skull_raw.nii.gz x -z y skull_FS.nii.gz
e mri_convert -odt uchar ./skull_FS.nii.gz ./skull_FS.nii.gz

# create skull surface from volume
e mri_tessellate ./skull_FS.nii.gz 255 ./skull_FS.fsmesh
# smooth skull surface
e mris_smooth -n 5 ./skull_FS.fsmesh ./skull_FS.fsmesh
# convert skull surface to stl format
e mris_convert ./skull_FS.fsmesh ./skull_FS.stl

# meshfix steps for skull surface
#e ${BINDIR}/meshfix skull_FS.stl -a 2.0 -u 5 --vertices $[$NUMBER_OF_VERTICES/3] -q -o skull.off
e ${BINDIR}/meshfix skull_FS.stl -a 2.0 -u 5 --vertices $NUMBER_OF_VERTICES -q -o skull.off
e ${BINDIR}/meshfix skull.off -a 2.0 -q -o skull.off # this steps included for final cleaning
e ${BINDIR}/meshfix skull.off -a 2.0 -u 1 -q -o skull.off
e ${BINDIR}/meshfix skull.off -a 2.0 -q -o skull.off

# decouple skull surface from csf (should not be necessary after volume decoupling; included just for safety)
while ${BINDIR}/meshfix skull.off $M2M_DIR/csf.stl --shells 2 --no-clean --intersect; do
 e ${BINDIR}/meshfix skull.off $M2M_DIR/csf.stl -a 2.0 --shells 2 --decouple-outout 0 -o skull.off
 e ${BINDIR}/meshfix skull.off $M2M_DIR/csf.stl -a 2.0 --shells 2 --cut-inner 0 -o skull.off
 e ${BINDIR}/meshfix skull.off -a 2.0 -u 1 -o skull.off
done

# final skull surface in stl format
e ${BINDIR}/meshfix skull.off --stl -q -o $M2M_DIR/skull.stl

# create volume mask from skull.stl
e ${BINDIR}/meshfix $M2M_DIR/skull.stl --no-clean --fsmesh -o skull.fsmesh
e mris_transform --dst ../ref_FS.nii.gz --src ../ref_FS.nii.gz ./skull.fsmesh ../unity.xfm ./skull.fsmesh
e mris_fill -r 0.5 -c ./skull.fsmesh $M2M_DIR/skull.nii.gz
e fslswapdim $M2M_DIR/skull.nii.gz x z -y $M2M_DIR/skull.nii.gz


# ==================================
SAY "Get skin mask from T1 and T2"
# ==================================
if $KEEP_MASKS; then
  SAY "using existing skin_raw volume!"
  if [ ! -f $MASK_DIR/skin_raw.nii.gz ]; then 
    echo "ERROR: Did not find $MASK_DIR/skin_raw.nii.gz"
    exit 1
  fi
  e imcp $MASK_DIR/skin_raw skin_raw

else
  if [ -f $M2M_DIR/T1.nii.gz ]; then
    e mv skull_raw_outskin_mask.nii.gz skin_raw.nii.gz
  else
    e mv csf_raw_outskin_mask.nii.gz skin_raw.nii.gz
  fi
  SKIN_IMAGE=$M2M_DIR/T1fs_nu_conform

  # get robust estimate of minimal skin intensity
  e fslmaths skin_raw -sub $M2M_DIR/skull -bin -mul $SKIN_IMAGE skin_estimate
  skin_threshold=`fslstats skin_estimate -P 10`
  e echo 'Threshold for skin' $skin_threshold

  # create enlarged skull volume
  e fslmaths $M2M_DIR/skull.nii.gz -dilM -bin skull_large.nii.gz

  # create better skin mask based on T1_conform (including volume decoupling from skull)
  e fslmaths $SKIN_IMAGE -s 1 -thr $skin_threshold -add skin_raw.nii.gz -add skull_large.nii.gz -bin skin_raw.nii.gz

  # get rid of regions outside of head (occuring for 1 ch Tx/Rx coil and T1-nu-image)
  e fslmaths skin_raw -ero -ero skin_raw_ero
  e cluster --in=skin_raw_ero --thresh=0.5 --no_table --osize=skin_raw_ero
  VAR=`fslstats skin_raw_ero -R`
  VAR=`echo $VAR | awk '{print $2}'`
  e fslmaths skin_raw_ero -thr $VAR -bin skin_raw_ero
  e fslmaths skin_raw_ero -dilM -dilM -dilM -dilM -mul skin_raw -bin -fillh skin_raw

  # test whether slightly smaller versions of skin_raw fit maximal intensity gradient in T1 better
  e fslmaths $M2M_DIR/T1fs_nu_conform -edge T1fs_nu_edge
  e imcp skin_raw skin_tmp_it1
  MAXEDGESTRTH=0
  BESTFITIT=1
  for (( i=1; i <= 5; i++ )); do
    e fslmaths skin_tmp_it$i -s 1 -thr 0.7 -bin skin_tmp_small
    e fslmaths skin_tmp_it$i -sub skin_tmp_small -mul T1fs_nu_edge skin_tmp_strength
    EDGESTRTH=`fslstats skin_tmp_strength -M`

    if (( $(bc <<< "$EDGESTRTH > $MAXEDGESTRTH") == 1 )); then
       MAXEDGESTRTH=$EDGESTRTH;
       BESTFITIT=$i;
    fi
    e immv skin_tmp_small 'skin_tmp_it'$((i+1))
  done
  e fslmaths 'skin_tmp_it'$BESTFITIT -add skull_large -bin skin_raw

  # store in mask_prep
  e imcp skin_raw $MASK_DIR/skin_raw
fi

# convert skin mask to FS space and format
e fslswapdim skin_raw.nii.gz x -z y skin_FS.nii.gz
e mri_convert -odt uchar ./skin_FS.nii.gz ./skin_FS.nii.gz

# create skin surface from volume
e mri_tessellate ./skin_FS.nii.gz 255 ./skin_FS.fsmesh
# smooth skin surface
e mris_smooth -n 5 ./skin_FS.fsmesh ./skin_FS.fsmesh
# convert skin surface to stl format
e mris_convert ./skin_FS.fsmesh ./skin_FS.stl

# meshfix steps for skin surface
e ${BINDIR}/meshfix skin_FS.stl -q --stl -o skin_FS.stl # work-around to prevent segmentation fault in next step that occured in one case
#e ${BINDIR}/meshfix skin_FS.stl -a 2.0 -u 5 --vertices $[$NUMBER_OF_VERTICES/4] -q -o skin.off
e ${BINDIR}/meshfix skin_FS.stl -a 2.0 -u 5 --vertices $NUMBER_OF_VERTICES -q -o skin.off
e ${BINDIR}/meshfix skin.off -a 2.0 -q -o skin.off # this steps included for final cleaning
e ${BINDIR}/meshfix skin.off -a 2.0 -u 1 -q -o skin.off
e ${BINDIR}/meshfix skin.off -a 2.0 -q -o skin.off

# decouple skin surface from skull (should not be necessary after volume decoupling; included just for safety)
while ${BINDIR}/meshfix skin.off $M2M_DIR/skull.stl --shells 2 --no-clean --intersect; do
 e ${BINDIR}/meshfix skin.off $M2M_DIR/skull.stl -a 2.0 --shells 2 --decouple-outout 0 -o skin.off
 e ${BINDIR}/meshfix skin.off $M2M_DIR/skull.stl -a 2.0 --shells 2 --cut-inner 0 -o skin.off
 e ${BINDIR}/meshfix skin.off -a 2.0 -u 1 -o skin.off
done

# final skin surface in stl format
e ${BINDIR}/meshfix skin.off --stl -q -o $M2M_DIR/skin.stl

# create volume mask from skin.stl
e ${BINDIR}/meshfix $M2M_DIR/skin.stl --no-clean --fsmesh -o skin.fsmesh

e mris_transform --dst ../ref_FS.nii.gz --src ../ref_FS.nii.gz ./skin.fsmesh ../unity.xfm ./skin.fsmesh
e mris_fill -r 0.5 -c ./skin.fsmesh $M2M_DIR/skin.nii.gz
e fslswapdim $M2M_DIR/skin.nii.gz x z -y $M2M_DIR/skin.nii.gz

e2 cd $OLD_PWD



####
# debugging (in m2m_<subjectID> directory):

# compare final wm and gm surfaces with original FreeSurfer output
# gmsh tmp/wm.lh.stl tmp/wm.rh.stl wm.stl &
# gmsh tmp/gm.lh.stl tmp/gm.rh.stl gm.stl &

# check final wm and gm surfaces
# gmsh wm.stl gm.stl &

# check registration to T1fs_conform
# fslview T1fs_conform T1_conform T2fs_conform T2_conform &

# check registration of T1fs_conform to MNI space
# fslview tmp/T1fs_conform_MNI.nii.gz & # load 1 mm MNI template

# check bet mask used as input to betsurf
# fslview T1fs_conform tmp/T1fs_conform_betted -l Blue-Lightblue&

# check final regions overlaid onto T1fs_conform
# fslview T1fs_conform  skin -l Red -t 0.4  skull -l Yellow -t 0.6  csf -l Blue -t 0.6  wm -l Yellow -t 0.7  gm -l Green -t 0.5  cerebellum -l Green -t 0.6  ventricles -l Blue -t 0.6 &

# check final stl surfaces
# gmsh wm.stl gm.stl csf.stl skull.stl skin.stl ventricles.stl cerebellum.stl &

