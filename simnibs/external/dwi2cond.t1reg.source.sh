
# check for mri2mesh or headreco output files
if [ ! -d $M2M_DIR ]; then e echo "ERROR: Directory $M2M_DIR doesn't exist or is not readable! Did you run mri2mesh or headreco before?"; exit; fi
for i in {$M2M_DIR/T1.nii.gz,$M2M_DIR/segmentation/labeling.nii.gz,$M2M_DIR/segmentation/T1_bias_corrected.nii.gz}; do
     if [ ! -f $i ]; then 
       e echo "ERROR: Input file '$i' doesn't exist. Did you run charm?"
       exit 1
     fi
done


# registration to structural T1
############################################
  e echo "registering results to T1fs_conform"
  OLD_PWD=`pwd`
  if [ ! -d $DTICONF_DIR ]; then e mkdir $DTICONF_DIR; fi
  e2 cd $DTICONF_DIR

# create T1 brain mask
  if [ ! -f $M2M_DIR'/final_tissues.nii.gz' ]; then e echo "ERROR: This version of dwi2cond requires a charm segmentation"; exit; fi
  e fslmaths $M2M_DIR/segmentation/labeling.nii.gz -thr 1 -uthr 499 -bin T1_brainmask
  e fslmaths $M2M_DIR/segmentation/T1_bias_corrected -mas T1_brainmask T1_brain
  e fslmaths T1_brainmask -edge -thr 0.3 -bin T1_brainrim_QA # for QA

# registration
  FLIRTDOF=12
  if [ $T1REGMTHD == '6dof' ]; then
	e echo 'using 6 DOF'
    FLIRTDOF=6
  fi
  e fslmaths $M2M_DIR/T1 -bin -s 1 T1_mask
  e fslcpgeom T1_brain T1_mask # work around:
  # in rare occasions, the dimensions of the two images seem not to be fully identical (header information
  # is the same when comparing manually; maybe some spurious differences in later digits?) so that flirt
  # throws an exception unless header information is copied
  e flirt -in $DTIRAW_DIR/DTI_FA -ref T1_brain -refweight T1_mask -omat FA2T1.mat -dof $FLIRTDOF
  if [ $T1REGMTHD == 'nonl' ]; then
    e echo 'using nonlinear registration'
    e fnirt --in=$DTIRAW_DIR/DTI_FA --ref=T1_brain --aff=FA2T1.mat --cout=FA2T1_warp --subsamp=8,4,2,2
    e vecreg -i $DTIRAW_DIR/DTI_tensor -o $M2M_DIR/DTI_coregT1_tensor -r T1_brain -w FA2T1_warp
  else
    e vecreg -i $DTIRAW_DIR/DTI_tensor -o $M2M_DIR/DTI_coregT1_tensor -r T1_brain -t FA2T1.mat
  fi
  e fslmaths $M2M_DIR/DTI_coregT1_tensor -mas T1_brainmask $M2M_DIR/DTI_coregT1_tensor
  

# preparing some results for QA
  e fslmaths $M2M_DIR/DTI_coregT1_tensor -tensor_decomp DTI_coregT1
  e rm DTI_coregT1.* DTI_coregT1_L?.* DTI_coregT1_M?.* DTI_coregT1_V2.* DTI_coregT1_V3.*

  e flirt -in $DTIRAW_DIR/DTI_FA -ref T1_brain -refweight T1_mask -omat FA2T1_QA.mat -out DTI_FA_6dof_QA -dof 6
  e fslmaths DTI_FA_6dof_QA -mas T1_brainmask DTI_FA_6dof_QA
  if [ -f $DTIRAW_DIR/DTI_sse.nii.gz ]; then 
    e flirt -in $DTIRAW_DIR/DTI_sse -ref T1_brain -applyxfm -init FA2T1_QA.mat -out DTI_SSE_6dof_QA
    e fslmaths DTI_SSE_6dof_QA -mas T1_brainmask DTI_SSE_6dof_QA
  fi

  if [ -d $RAW_DIR ]; then
    e flirt -in $RAW_DIR/DWIraw_FA -ref T1_brain -refweight T1_mask -omat FA2T1_QA.mat -out DTIraw_FA_6dof_QA -dof 6
    e flirt -in $RAW_DIR/DWIraw_SSE -ref T1_brain -applyxfm -init FA2T1_QA.mat -out DTIraw_SSE_6dof_QA
    e fslmaths DTIraw_FA_6dof_QA -mas T1_brainmask DTIraw_FA_6dof_QA
    e fslmaths DTIraw_SSE_6dof_QA -mas T1_brainmask DTIraw_SSE_6dof_QA
  fi

  e nii2msh -ev $M2M_DIR/DTI_coregT1_tensor.nii.gz $M2M_DIR/${SUBJECT}.msh $D2C_DIR/first_ev_for_check.msh
  e rm T1_mask.* T1_brainmask.* FA2T1_QA.mat
  if [ $TIDY_UP ]; then
    for i in {$RAW_DIR,$TOPUP_DIR,$EDDY_DIR,$FM_DIR,$EDDYCOR_DIR}; do
        if [ -d $i ]; then
            e rm $i/*
            e rmdir $i
        fi
    done
  fi
  e2 cd $OLD_PWD

