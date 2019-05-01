# check for first input file
if [ ! -d $M2M_DIR -o ! -f $T1fs ]; then echo "ERROR: Input file doesn't exist. Did you call this script from mri2mesh?"; exit 1; fi;

# test: do we already have the recon-all results?
for i in $FS_DIR/{surf/{l,r}h.{white,pial},mri/{aparc.a2009s+aseg,orig,rawavg}.mgz}; do if [ ! -f $i ]; then run_brainf=true; fi; done;

if $run_brainf; then
	# ==================================
	SAY "running FreeSurfer recon-all on T1fs!"
	# ==================================

	# ensure standard orientation of all input images
	e fslmaths $T1fs $M2M_DIR/T1fs # this is needed as fslreorient2std ignores FSLOUTPUTTYPE
	e fslreorient2std $M2M_DIR/T1fs $M2M_DIR/T1fs
	if [ -n "$T1" ]; then e fslmaths $T1 $M2M_DIR/T1; e fslreorient2std $M2M_DIR/T1 $M2M_DIR/T1; fi
	if [ -n "$T2" ]; then e fslmaths $T2 $M2M_DIR/T2; e fslreorient2std $M2M_DIR/T2 $M2M_DIR/T2; fi
	if [ -n "$T2fs" ]; then e fslmaths $T2fs $M2M_DIR/T2fs; e fslreorient2std $M2M_DIR/T2fs $M2M_DIR/T2fs; fi
	
	# ensure radiological orientation of all input images
	for i in T1 T1fs T2 T2fs; do
	  # ensuring RADIOLOGICAL orientation
	  if [ ! -f $M2M_DIR/$i.nii.gz ]; then continue; fi;
	  VAR=`fslorient -getorient $M2M_DIR/$i.nii.gz`;
	  if [ $VAR = NEUROLOGICAL ] ; then
	      	e echo $i'.nii.gz: enforcing radiological orientation'
	      	e fslswapdim $M2M_DIR/$i.nii.gz -x y z $M2M_DIR/$i.nii.gz
	      	e fslorient -forceradiological $M2M_DIR/$i.nii.gz
	  fi
	done;

	# preprocess T1fs
	# resample to desired voxel size
	VOXSIZE=1 # for now, freesurfer 5.3.0 works only properly with 1 mm iso voxel

	e fslswapdim $M2M_DIR/T1fs x -z y $M2M_DIR/tmpImg # get to LIA, otherwise mri_convert will cause a shift in already conformed images
        e mri_convert -cs $VOXSIZE $M2M_DIR/tmpImg.nii.gz $M2M_DIR/tmpImg.nii.gz
        e mri_convert -cm $M2M_DIR/tmpImg.nii.gz $M2M_DIR/tmpImg.nii.gz # work-around to get the final image size for vox sizes < 1 mm
	df_conf=`fslval $M2M_DIR/tmpImg dim1`
        e rm $M2M_DIR/tmpImg.nii.gz

	e fslswapdim $M2M_DIR/T1fs x -z y $M2M_DIR/tmpImg
	e mri_convert -nc -cs $VOXSIZE -oni $df_conf -onj $df_conf -onk $df_conf -rt cubic $M2M_DIR/tmpImg.nii.gz $M2M_DIR/T1fs_resamp.nii.gz
        # -nc flag ensures that intensities are kept (for t1t2ratio)
	e fslswapdim $M2M_DIR/T1fs_resamp x z -y $M2M_DIR/T1fs_resamp # get to LAS again
	e rm $M2M_DIR/tmpImg.nii.gz

	# set header so that qform and sform are identical to voxel-to-surface mapping used by freesurfer and simnibs
        offset_x=`CALC $VOXSIZE*$df_conf/2`
	offset_y=`CALC -$VOXSIZE*$df_conf/2`
	offset_z=`CALC -$VOXSIZE*$df_conf/2+$VOXSIZE`
        e fslorient -setqform -$VOXSIZE 0 0 $offset_x 0 $VOXSIZE 0 $offset_y 0 0 $VOXSIZE $offset_z 0 0 0 1 $M2M_DIR/T1fs_resamp
        e fslorient -setsform -$VOXSIZE 0 0 $offset_x 0 $VOXSIZE 0 $offset_y 0 0 $VOXSIZE $offset_z 0 0 0 1 $M2M_DIR/T1fs_resamp

	# get head size in z-direction (registered MNI is 0 outside head, so that fslstats works)
	e flirt -in $STD_MNI -ref $M2M_DIR/T1fs_resamp -out $M2M_DIR/tmp/MNIfs_reg
	var=`fslstats $M2M_DIR/tmp/MNIfs_reg -w`
	zmin=`echo $var | awk '{print $5}'`
	zsize=$(${BINDIR}/expr $df_conf - $zmin)	

	# cut neck, pad zeros to get back to original dimensions, and fix header
	e fslroi $M2M_DIR/T1fs_resamp $M2M_DIR/T1fs_roi 0 -1 0 -1 $zmin $zsize
	e fslcreatehd $df_conf $df_conf $zmin 1 $VOXSIZE $VOXSIZE $VOXSIZE 1 0 0 0 16 $M2M_DIR/tmpImg
	e fslmerge -z $M2M_DIR/T1fs_roi $M2M_DIR/tmpImg $M2M_DIR/T1fs_roi
	e rm $M2M_DIR/tmpImg.nii.gz
        e fslcpgeom $M2M_DIR/T1fs_resamp $M2M_DIR/T1fs_roi

	# input to freesurfer has to be LIA, otherwise mri_convert inside freesurfer will change qform and sform
	e fslswapdim $M2M_DIR/T1fs_roi x -z y $M2M_DIR/tmp/T1fs_roi_FS
     
	# Start surface reconstuction and segmentation in freesurfer
	# (this creates the folder structure $SUBJECT/{mri,surf,...} )
	if [ $T2MASK = 1 ]; then
	  e echo 'using T2 to mask out dura'
	  e fslswapdim $M2M_DIR/T2 x -z y $M2M_DIR/tmp/T2_FS
	  e flirt -in $M2M_DIR/tmp/T2_FS -ref $M2M_DIR/tmp/T1fs_roi_FS -interp spline -cost mutualinfo -searchcost mutualinfo \
	  -dof 6 -out $M2M_DIR/tmp/T2_regT1fs_roi

	  VAR=`fslstats $M2M_DIR/tmp/T2_regT1fs_roi -M`
	  e fslmaths $M2M_DIR/tmp/T2_regT1fs_roi -thr $VAR -bin $M2M_DIR/tmp/T2_regT1fs_roi_mask

	  e cluster --in=$M2M_DIR/tmp/T2_regT1fs_roi_mask --thresh=0.5 --no_table --osize=$M2M_DIR/tmp/T2_regT1fs_roi_mask
	  VAR=`fslstats $M2M_DIR/tmp/T2_regT1fs_roi_mask -R`
	  VAR=`echo $VAR | awk '{print $2}'`
	  e fslmaths $M2M_DIR/tmp/T2_regT1fs_roi_mask -thr $VAR -bin -fillh $M2M_DIR/tmp/T2_regT1fs_roi_mask

	  e fslmaths $M2M_DIR/tmp/T1fs_roi_FS -mul $M2M_DIR/tmp/T2_regT1fs_roi_mask -add $M2M_DIR/tmp/T1fs_roi_FS \
	  -div 2 $M2M_DIR/T1fs_roi_T2masked_FS

	  e recon-all -i $M2M_DIR/T1fs_roi_T2masked_FS.nii.gz -s $FS_DIR
	  T2PIAL=0
	else
	  e recon-all -i $M2M_DIR/tmp/T1fs_roi_FS.nii.gz -s $FS_DIR
	fi

	if [ $T2PIAL = 0 ]; then
 	  e recon-all -s $FS_DIR -all -3T -cubic -cm
	else
	  e echo 'using T2 to improve segmentation of pial surface'

	  # preprocess T2
	  # get head size in z-direction (registered MNI is 0 outside head, so that fslstats works)
	  e flirt -in $STD_MNI -ref $M2M_DIR/T2.nii.gz -out $M2M_DIR/tmp/MNIfs_regT2 
	  var=`fslstats $M2M_DIR/tmp/MNIfs_regT2 -w`
	  zmin=`echo $var | awk '{print $5}'`
	  zsize=`echo $var | awk '{print $6}'`
	  # get dimensions of T2
	  df1=`fslval $M2M_DIR/T2.nii.gz dim1`
	  df2=`fslval $M2M_DIR/T2.nii.gz dim2`
	  # cut neck
	  e fslroi $M2M_DIR/T2.nii.gz $M2M_DIR/T2_roi 0 $df1 0 $df2 $zmin $zsize

	  e recon-all -T2 $M2M_DIR/T2_roi.nii.gz -s $FS_DIR
	  e recon-all -s $FS_DIR -all -3T -cubic -T2pial -cm
	fi;
else SAY "skipping recon-all, using previous results to save time!"; fi


# ==================================
SAY "Convert WM and GM hemispheres to STL"
# ==================================
# Convert the surface meshes of the hemispheres to STL
for j in {r,l}h; do
  e mris_convert $FS_DIR/surf/$j.white $M2M_DIR/tmp/wm.$j.stl
  e mris_convert $FS_DIR/surf/$j.pial  $M2M_DIR/tmp/gm.$j.stl
done

# delete links created by freesurfer
for i in {fsaverage,lh.EC_average,rh.EC_average}; do
  if [ -L $i ]; then rm $i; fi
done

OLD_PWD=`pwd`
e2 cd $M2M_DIR

# get filter kernel for broadening the corpus callosum in x-direction
e cp $TEMPLATEPATH/CC_filterkernel.nii.gz $M2M_DIR/tmp/.

# transform T1fs and T1fs_nu from conformed FS space to FSL
e mri_convert $FS_DIR/mri/orig.mgz $M2M_DIR/T1fs_conform.nii.gz
e fslswapdim $M2M_DIR/T1fs_conform.nii.gz x z -y $M2M_DIR/T1fs_conform.nii.gz
e cp $M2M_DIR/T1fs_conform.nii.gz $SUBJECTS_DIR/$SUBJECT"_T1fs_conform.nii.gz"

e mri_convert $FS_DIR/mri/nu.mgz $M2M_DIR/T1fs_nu_conform.nii.gz
e fslswapdim $M2M_DIR/T1fs_nu_conform.nii.gz x z -y $M2M_DIR/T1fs_nu_conform.nii.gz


# ===========================================================
SAY "creating subcortical mask for joining the hemispheres"
# ===========================================================
LABELS_SUBCORT=(4 10 14 24 28 31 43 49 60 63)

e2 cd $M2M_DIR/tmp

# convert FS label mask to FSL
e mri_convert $FS_DIR/mri/aparc.a2009s+aseg.mgz ./segmentation-mask.nii.gz
e fslswapdim segmentation-mask.nii.gz x z -y segmentation-mask.nii.gz

# create CC mask 
e fslmaths segmentation-mask.nii.gz -thr 251 -uthr 255 -bin CC.nii.gz
e fslmaths CC.nii.gz -kernel file CC_filterkernel.nii.gz -fmean CC_filt.nii.gz -odt float
e fslmaths CC_filt.nii.gz -thr 0.1 -bin CC_final.nii.gz

# create mask with shifted central CSF slab (used to mask the fornix)
e2 cat > yshift.mat <<EOF
     1     0     0     0
     0     1     0    10
     0     0     1     0
     0     0     0     1
EOF

e fslmaths segmentation-mask.nii.gz -thr 24 -uthr 24 -dilM -bin CSFslab.nii.gz
e flirt -in CSFslab.nii.gz -ref CSFslab.nii.gz -applyxfm -init yshift.mat -out CSFslab_shift.nii.gz

# create large subcortical mask (used to join hemispheres again)
e fslmaths segmentation-mask.nii.gz -mul 0 -bin subcortical.nii.gz
for i in ${LABELS_SUBCORT[*]}; do
 e fslmaths segmentation-mask.nii.gz -thr $i -uthr $i -add subcortical.nii.gz -bin subcortical.nii.gz
done
e fslmaths subcortical.nii.gz -add CC_filt.nii.gz -add CSFslab_shift.nii.gz -dilM -ero -bin subcortical.nii.gz

# convert subcortical mask to FS space
e fslswapdim subcortical.nii.gz x -z y subcortical_FS.nii.gz 
e mri_convert -odt uchar ./subcortical_FS.nii.gz ./subcortical_FS.nii.gz

# store "reference" volume with header information for later usage with mris_transform
e cp subcortical_FS.nii.gz ../ref_FS.nii.gz

# create subcortical surface from volume, smooth and convert to stl
e mri_tessellate ./subcortical_FS.nii.gz 255 ./subcortical_FS.fsmesh
e mris_smooth -n 5 ./subcortical_FS.fsmesh ./subcortical_FS.fsmesh
e mris_convert ./subcortical_FS.fsmesh ./subcortical_FS.stl

# create optic radiation mask and surface, smooth surface and convert to stl
e fslmaths segmentation-mask.nii.gz -thr 85 -uthr 85 -bin opticrad.nii.gz

e fslswapdim opticrad.nii.gz x -z y opticrad_FS.nii.gz 
e mri_convert -odt uchar ./opticrad_FS.nii.gz ./opticrad_FS.nii.gz

e mri_tessellate ./opticrad_FS.nii.gz 255 ./opticrad_FS.fsmesh
e mris_smooth -n 5 ./opticrad_FS.fsmesh ./opticrad_FS.fsmesh
e mris_convert ./opticrad_FS.fsmesh ./opticrad_FS.stl


#############################
# meshfix steps start below
#############################

# clean subcortical surface
e ${BINDIR}/meshfix subcortical_FS.stl -a 2.0 --remove-handles -q -o subcortical.off
e ${BINDIR}/meshfix subcortical.off -a 2.0 -u 5 --vertices $[$NUMBER_OF_VERTICES/8] -q -o subcortical.off
# create large, small and very small subcortical surfaces
e ${BINDIR}/meshfix subcortical.off -a 2.0 --dilate -1 -q -o subcortical_small
e ${BINDIR}/meshfix subcortical.off -a 2.0 --dilate -1 -q -o subcortical_small_small
e ${BINDIR}/meshfix subcortical.off -a 2.0 --dilate  1 -q -o subcortical_large
e ${BINDIR}/meshfix subcortical_small.off -a 2.0 -u 5 -q -o subcortical_small

# clean optic radition surface
e ${BINDIR}/meshfix opticrad_FS.stl -a 2.0 --remove-handles -q -o opticrad.off
e ${BINDIR}/meshfix opticrad.off -a 2.0 -u 5 -q -o opticrad.off
e ${BINDIR}/meshfix opticrad.off -a 2.0 --dilate  1 -q -o opticrad_large
e ${BINDIR}/meshfix opticrad_large.off -a 2.0 --dilate  1 -q -o opticrad_large


# ======================================
SAY "preparing wm and gm hemispheres"
# ======================================
for i in {r,l}h; do
 e ${BINDIR}/meshfix wm.$i.stl -a 2.0 --remove-handles -q -o wm.${i}_fixed
 e ${BINDIR}/meshfix wm.${i}_fixed.off -a 2.0 -u 5 --vertices $NUMBER_OF_VERTICES -q -o wm.${i}_fixed

 e ${BINDIR}/meshfix wm.${i}_fixed.off -a 2.0 -u 1 -q -o wm.${i}_fixed
 e ${BINDIR}/meshfix wm.${i}_fixed.off -a 2.0 -q -o wm.${i}_fixed

 e ${BINDIR}/meshfix gm.$i.stl -a 2.0 -u 5 --vertices $[115*$NUMBER_OF_VERTICES/100] --remove-handles -q -o gm.${i}_fixed
 e ${BINDIR}/meshfix gm.${i}_fixed.off -a 2.0 -u 5 -q -o gm.${i}_fixed

 e ${BINDIR}/meshfix gm.${i}_fixed.off -a 2.0 -u 1 -q -o gm.${i}_fixed
 e ${BINDIR}/meshfix gm.${i}_fixed.off -a 2.0 -q -o gm.${i}_fixed
done

# cutting optic radiations away
for i in {r,l}h; do
  e ${BINDIR}/meshfix wm.${i}_fixed.off opticrad_large.off -a 2.0 --shells 2 --cut-inner 0 -q -o wm.${i}_fixed.off
  e ${BINDIR}/meshfix gm.${i}_fixed.off opticrad_large.off -a 2.0 --shells 2 --cut-inner 0 -q -o gm.${i}_fixed.off
done


# ======================================
SAY "joining wm hemispheres"
# ======================================
# join wm rh with small subcortical surface
e ${BINDIR}/meshfix wm.rh_fixed.off subcortical_small.off -a 2.0 --shells 2 -j -o wm.rh_fixed2
# join wm rh with wm lh
e ${BINDIR}/meshfix wm.rh_fixed2.off wm.lh_fixed.off -a 2.0 --shells 2 -j -o wm_fixed

# WM: remove throats that can occur at the joining lines
e cp wm_fixed.off wm_fixed2.off
WMONE=0
for i in {1..5}; do
 if ${BINDIR}/meshfix wm_fixed2.off subcortical_small.off --shells 2 --no-clean --intersect; then
    e ${BINDIR}/meshfix wm_fixed2.off subcortical_small.off -a 2.0 --shells 2 --decouple-outout 0 -o wm_fixed2.off
    e ${BINDIR}/meshfix wm_fixed2.off subcortical_small.off -a 2.0 --shells 2 --cut-inner 0 -o wm_fixed2.off
    let WMONE++;
 else break;
 fi
done

# WM: push outside so that WM is as large as original WM
WMTWO=0
for i in {1..5}; do
 if ${BINDIR}/meshfix wm_fixed2.off subcortical.off --shells 2 --no-clean --intersect; then
    e ${BINDIR}/meshfix wm_fixed2.off subcortical.off -a 2.0 --shells 2 --decouple-outout 0 -o wm_fixed2.off
    let WMTWO++;
 else break;
 fi
done


# ======================================
SAY "joining gm hemispheres"
# ======================================
# join gm rh with subcortical surface
e ${BINDIR}/meshfix gm.rh_fixed.off subcortical_large.off -a 2.0 --shells 2 -j -o gm.rh_fixed2
# join gm rh with gm lh 
e ${BINDIR}/meshfix gm.rh_fixed2.off gm.lh_fixed.off -a 2.0 --shells 2 -j -o gm_fixed

# GM: remove throats that can occur at the joining lines
e cp gm_fixed.off gm_fixed2.off
GMONE=0
for i in {1..5}; do
 if ${BINDIR}/meshfix gm_fixed2.off subcortical.off --shells 2 --no-clean --intersect; then
    e ${BINDIR}/meshfix gm_fixed2.off subcortical.off -a 2.0 --shells 2 --decouple-outout 0 -o gm_fixed2.off
    e ${BINDIR}/meshfix gm_fixed2.off subcortical.off -a 2.0 --shells 2 --cut-inner 0 -o gm_fixed2.off
    let GMONE++;
 else break;
 fi
done

# GM: push outside so that GM is larger than WM
GMTWO=0
for i in {1..5}; do
 if ${BINDIR}/meshfix gm_fixed2.off subcortical_large.off --shells 2 --no-clean --intersect; then
    e ${BINDIR}/meshfix gm_fixed2.off subcortical_large.off -a 2.0 --shells 2 --decouple-outout 0 -o gm_fixed2.off
    let GMTWO++;
 else break;
 fi
done


# ======================================
SAY "decoupling final gm from final wm"
# ======================================
# initial cleaning of surfaces before decoupling
e ${BINDIR}/meshfix gm_fixed2.off -a 2.0 -u 1 -q -o gm_fixed3.off
e ${BINDIR}/meshfix gm_fixed3.off -a 2.0 -q -o gm_fixed3.off
e ${BINDIR}/meshfix gm_fixed3.off -a 2.0 -u 1 -q -o gm_fixed3.off
e ${BINDIR}/meshfix gm_fixed3.off -a 2.0 -q -o gm_fixed3.off

e ${BINDIR}/meshfix wm_fixed2.off -a 2.0 -u 1 -q -o wm_fixed3.off
e ${BINDIR}/meshfix wm_fixed3.off -a 2.0 -q -o wm_fixed3.off
e ${BINDIR}/meshfix wm_fixed3.off -a 2.0 -u 1 -q -o wm_fixed3.off
e ${BINDIR}/meshfix wm_fixed3.off -a 2.0 -q -o wm_fixed3.off

# intial decoupling for two iterations (push GM out)
K=0
for i in {1..2}; do
 if ${BINDIR}/meshfix gm_fixed3.off wm_fixed3.off --shells 2 --no-clean --intersect; then
  let K++;
  e echo '** GM-WM decoupling: Iteration '$K' **'

  # push GM out (3x), cut parts within WM if still necessary
  for i in {1..3}; do
    if ${BINDIR}/meshfix gm_fixed3.off wm_fixed3.off --shells 2 --no-clean --intersect; then
       e ${BINDIR}/meshfix gm_fixed3.off wm_fixed3.off -a 2.0 --shells 2 --decouple-outout 0 -o gm_fixed3.off
    else break;
    fi
  done
  e ${BINDIR}/meshfix gm_fixed3.off wm_fixed3.off -a 2.0 --shells 2 --cut-inner 0 -o gm_fixed3.off

  # cleaning of gm after decoupling steps
  e ${BINDIR}/meshfix gm_fixed3.off -a 2.0 -u 1 -q -o gm_fixed3.off
  e ${BINDIR}/meshfix gm_fixed3.off -a 2.0 -q -o gm_fixed3.off
 else break;
 fi
done

# after decoupling most of the intersections, fine-tune for volume meshing
e ${BINDIR}/meshfix gm_fixed3.off gm_fixed3.off --fineTuneOut 0.2 4 --shells 2 -o gm_fixed3.off
e ${BINDIR}/meshfix gm_fixed3.off -a 2.0 -u 1 -q -o gm_fixed3.off
e ${BINDIR}/meshfix gm_fixed3.off -a 2.0 -q -o gm_fixed3.off

e ${BINDIR}/meshfix wm_fixed3.off gm_fixed3.off --fineTuneIn  0.2 4 --shells 2 -o wm_fixed3.off
e ${BINDIR}/meshfix wm_fixed3.off -a 2.0 -u 1 -q -o wm_fixed3.off
e ${BINDIR}/meshfix wm_fixed3.off -a 2.0 -q -o wm_fixed3.off

# remove WM spikes if necessary
if ${BINDIR}/meshfix gm_fixed3.off wm_fixed3.off --shells 2 --no-clean --intersect; then
  e echo 'cutting WM and decoupling it from GM to get rid of spurious WM spikes'

  # cut WM spikes
  e ${BINDIR}/meshfix wm_fixed3.off gm_fixed3.off -a 2.0 --shells 2 --cut-outer 0 -o wm_fixed3.off
  e ${BINDIR}/meshfix wm_fixed3.off -a 2.0 -u 1 -q -o wm_fixed3.off
  e ${BINDIR}/meshfix wm_fixed3.off -a 2.0 -q -o wm_fixed3.off

  # push WM spikes inside
  for i in {1..3}; do
   if ${BINDIR}/meshfix gm_fixed3.off wm_fixed3.off --shells 2 --no-clean --intersect; then
      e ${BINDIR}/meshfix wm_fixed3.off gm_fixed3.off -a 2.0 --shells 2 --decouple-inin 0 -o wm_fixed3.off
      # cleaning of wm after decoupling steps
      e ${BINDIR}/meshfix wm_fixed3.off -a 2.0 -u 1 -q -o wm_fixed3.off
      e ${BINDIR}/meshfix wm_fixed3.off -a 2.0 -q -o wm_fixed3.off
   else break;
   fi
  done
fi

# final decoupling (push GM out)
while ${BINDIR}/meshfix gm_fixed3.off wm_fixed3.off --shells 2 --no-clean --intersect; do
 let K++;
 e echo '** GM-WM decoupling: Iteration '$K' **'
 
 for i in {1..3}; do
  if ${BINDIR}/meshfix gm_fixed3.off wm_fixed3.off --shells 2 --no-clean --intersect; then
     e ${BINDIR}/meshfix gm_fixed3.off wm_fixed3.off -a 2.0 --shells 2 --decouple-outout 0 -o gm_fixed3.off
  else break;
  fi
 done
 e ${BINDIR}/meshfix gm_fixed3.off wm_fixed3.off -a 2.0 --shells 2 --cut-inner 0 -o gm_fixed3.off

 # cleaning of gm after decoupling steps
 e ${BINDIR}/meshfix gm_fixed3.off -a 2.0 -u 1 -q -o gm_fixed3.off
 e ${BINDIR}/meshfix gm_fixed3.off -a 2.0 -q -o gm_fixed3.off
done
e echo 'Overall iterations needed for preparing WM and GM:'$WMONE' '$WMTWO' '$GMONE' '$GMTWO' '$K

e ${BINDIR}/meshfix gm_fixed3.off --stl -q -o $M2M_DIR/gm.stl
e ${BINDIR}/meshfix wm_fixed3.off --stl -q -o $M2M_DIR/wm.stl


# ===========================================================
SAY "creating volume masks from final GM and WM surfaces"
# ===========================================================

# converting from STL to FreeSurfer mesh (using ${BINDIR}/meshfix, as mris_convert is VERY slow)
e ${BINDIR}/meshfix $M2M_DIR/wm.stl --no-clean --fsmesh -o wm.fsmesh
e ${BINDIR}/meshfix $M2M_DIR/gm.stl --no-clean --fsmesh -o gm.fsmesh

# unity transformation (FS format)
e2 cat > unity.xfm <<EOF
MNI Transform File
% tkregister2

Transform_Type = Linear;
Linear_Transform =
   1.0    0.0    0.0  0.0 
   0.0    1.0    0.0  0.0 
   0.0    0.0    1.0  0.0 ;
EOF

e cp unity.xfm ../unity.xfm

# apply unity trafos in conform space, as mris_fill cuts the volume masks otherwise (FS bug?)
e mris_transform --dst ./subcortical_FS.nii.gz --src ./subcortical_FS.nii.gz ./wm.fsmesh ./unity.xfm ./wm.fsmesh
e mris_transform --dst ./subcortical_FS.nii.gz --src ./subcortical_FS.nii.gz ./gm.fsmesh ./unity.xfm ./gm.fsmesh

# create volume masks
e mris_fill -c ./wm.fsmesh $M2M_DIR/wm.nii.gz
e mris_fill -c ./gm.fsmesh $M2M_DIR/gm.nii.gz

# transform to FSL space
e fslswapdim $M2M_DIR/wm.nii.gz x z -y $M2M_DIR/wm.nii.gz
e fslswapdim $M2M_DIR/gm.nii.gz x z -y $M2M_DIR/gm.nii.gz

e2 cd $OLD_PWD
