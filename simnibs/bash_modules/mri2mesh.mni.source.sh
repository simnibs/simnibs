# check for first input file
if [ ! -d $M2M_DIR -o ! -f $M2M_DIR/T1fs_nu_conform.nii.gz ]; then echo "ERROR: Input file doesn't exist. Did you call this script from mri2mesh?"; exit 1; fi;
if [ ! -d $M2M_DIR/toMNI ]; then mkdir $M2M_DIR/toMNI; fi


#############################
# create linear mni to conformed space trafos for electrode placement
#############################
e echo 'linear registrations to MNI (fsl flirt) ... '
e flirt -dof 12 -in $STD_MNI -ref $M2M_DIR/T1fs_nu_conform -out $M2M_DIR/tmp/MNI2conform_12DOF -omat $M2M_DIR/tmp/MNI2conform_12DOF_fsl.mat
e flirt -dof 6  -in $STD_MNI -ref $M2M_DIR/T1fs_nu_conform -out $M2M_DIR/tmp/MNI2conform_6DOF  -omat $M2M_DIR/tmp/MNI2conform_6DOF_fsl.mat

# convert matrices from FSL format to "SPM" format (NOTE: this works only correctly as long as both images are in LAS space)
# the latter trafos from conformed mm world space to MNI mm world space (and inverse)
VAR=`fslorient -getorient $STD_MNI`;
if [ $VAR = NEUROLOGICAL ] ; then echo e "ERROR: MNI template has to be in radiological convention!"; exit 1; fi;


m=(`fslorient -getqform $M2M_DIR/T1fs_nu_conform`)
e2 cat > $M2M_DIR/tmp/qform_conf.mat <<EOF
  ${m[0]} ${m[1]} ${m[2]} ${m[3]}
  ${m[4]} ${m[5]} ${m[6]} ${m[7]}
  ${m[8]} ${m[9]} ${m[10]} ${m[11]}
  ${m[12]} ${m[13]} ${m[14]} ${m[15]}
EOF

m=(`fslorient -getqform $STD_MNI`)
e2 cat > $M2M_DIR/tmp/qform_mni.mat <<EOF
  ${m[0]} ${m[1]} ${m[2]} ${m[3]}
  ${m[4]} ${m[5]} ${m[6]} ${m[7]}
  ${m[8]} ${m[9]} ${m[10]} ${m[11]}
  ${m[12]} ${m[13]} ${m[14]} ${m[15]}
EOF

# create internal FSL matrix for MNI template
VOXSIZE_X=`fslval $STD_MNI pixdim1`
VOXSIZE_Y=`fslval $STD_MNI pixdim2`
VOXSIZE_Z=`fslval $STD_MNI pixdim3`
e2 cat > $M2M_DIR/tmp/fsl_mni.mat <<EOF
  $VOXSIZE_X     0           0        0
     0        $VOXSIZE_Y     0        0
     0           0       $VOXSIZE_Z   0
     0           0           0        1
EOF

# create internal FSL matrix for conformed image
VOXSIZE_X=`fslval $M2M_DIR/T1fs_nu_conform pixdim1`
VOXSIZE_Y=`fslval $M2M_DIR/T1fs_nu_conform pixdim2`
VOXSIZE_Z=`fslval $M2M_DIR/T1fs_nu_conform pixdim3`
e2 cat > $M2M_DIR/tmp/fsl_conf.mat <<EOF
  $VOXSIZE_X        0          0        0
     0          $VOXSIZE_Y     0        0
     0              0      $VOXSIZE_Z   0
     0              0          0        1
EOF

# mm-to-mm mapping from MNI to conformed space (linear):
# Pos_conf = qform_conf*inv(fsl_conf)*R_fsl*fsl_mni*inv(qform_mni)*Pos_MNI
e convert_xfm -omat $M2M_DIR/tmp/qform_mni_inv.mat -inverse $M2M_DIR/tmp/qform_mni.mat
e convert_xfm -omat $M2M_DIR/tmp/postmat.mat -concat $M2M_DIR/tmp/fsl_mni.mat $M2M_DIR/tmp/qform_mni_inv.mat

e convert_xfm -omat $M2M_DIR/tmp/fsl_conf_inv.mat -inverse $M2M_DIR/tmp/fsl_conf.mat
e convert_xfm -omat $M2M_DIR/tmp/premat.mat -concat $M2M_DIR/tmp/qform_conf.mat $M2M_DIR/tmp/fsl_conf_inv.mat

e convert_xfm -omat $M2M_DIR/toMNI/MNI2conform_12DOF.mat -concat $M2M_DIR/tmp/MNI2conform_12DOF_fsl.mat $M2M_DIR/tmp/postmat.mat
e convert_xfm -omat $M2M_DIR/toMNI/MNI2conform_12DOF.mat -concat $M2M_DIR/tmp/premat.mat $M2M_DIR/toMNI/MNI2conform_12DOF.mat

e convert_xfm -omat $M2M_DIR/toMNI/MNI2conform_6DOF.mat -concat $M2M_DIR/tmp/MNI2conform_6DOF_fsl.mat $M2M_DIR/tmp/postmat.mat
e convert_xfm -omat $M2M_DIR/toMNI/MNI2conform_6DOF.mat -concat $M2M_DIR/tmp/premat.mat $M2M_DIR/toMNI/MNI2conform_6DOF.mat



#############################
# create nonlinear conformed to MNI space trafos
#############################
e echo 'non-linear registrations to MNI (fsl fnirt) ... takes some time'
e convert_xfm -omat $M2M_DIR/tmp/Conform2MNI_12DOF_fsl.mat -inverse $M2M_DIR/tmp/MNI2conform_12DOF_fsl.mat

e imcp $M2M_DIR/T1fs_nu_conform $M2M_DIR/tmp/T1fs_for_fnirt
if [ $MNIMASKSKULL = 1 ]; then
  e echo 'downweighting skull region for FNIRT registration'
  e fslmaths $M2M_DIR/gm -add $M2M_DIR/cerebellum -dilM -dilM -fillh -ero $M2M_DIR/tmp/brain_mask
  # note: gm mask from FreeSurfer is more robust than csf mask, which tends to be also incorrect for problematic T1w images with bright spongy bone
  e fslmaths $M2M_DIR/skin -sub $M2M_DIR/skull -add $M2M_DIR/tmp/brain_mask -binv -mul 10 $M2M_DIR/tmp/skull_masked
  e fslmaths $M2M_DIR/skin -sub $M2M_DIR/skull -add $M2M_DIR/tmp/brain_mask -bin -mul $M2M_DIR/tmp/T1fs_for_fnirt -add $M2M_DIR/tmp/skull_masked $M2M_DIR/tmp/T1fs_for_fnirt
fi
e fnirt --in=$M2M_DIR/tmp/T1fs_for_fnirt --aff=$M2M_DIR/tmp/Conform2MNI_12DOF_fsl.mat --config=T1_2_MNI152_2mm --refmask=$TEMPLATEPATH/MNI152_T1_2mm_head_dil --fout=$M2M_DIR/tmp/Conform2MNI_nonl_fsl_field


# getting transformation from MNI mm world space to conformed mm world space
e convert_xfm -omat $M2M_DIR/tmp/premat_inv.mat -inverse $M2M_DIR/tmp/premat.mat
e convertwarp --ref=$M2M_DIR/tmp/Conform2MNI_nonl_fsl_field --premat=$M2M_DIR/tmp/premat_inv.mat --warp1=$M2M_DIR/tmp/Conform2MNI_nonl_fsl_field --absout --out=$M2M_DIR/toMNI/MNI2Conform_nonl
# note: fsl stores the inverse warp (i.e. Conform2MNI_nonl_fsl_field stores the 
# warp from MNI-to-conform); the field coordinates in Conform2MNI_nonl_fsl_field therefore
# correspond to the transformation from MNI to "conformed FSL mm space"; concatenating the inverse trafo
# from conformed mm world space to conformed FSL mm space as "premat" (effectively as postmat in the normal
# transformation direction) thus results in a trafo from MNI mm world space to conformed mm world space


# getting transformation from conformed mm world space to MNI mm world space
e fslmaths $M2M_DIR/tmp/T1fs_for_fnirt -subsamp2 $M2M_DIR/tmp/T1fs_for_fnirt_ds2
e invwarp -w $M2M_DIR/tmp/Conform2MNI_nonl_fsl_field -o $M2M_DIR/tmp/MNI2Conform_nonl_fsl_field -r $M2M_DIR/tmp/T1fs_for_fnirt_ds2 
e convertwarp --ref=$M2M_DIR/tmp/MNI2Conform_nonl_fsl_field --premat=$M2M_DIR/tmp/postmat.mat --warp1=$M2M_DIR/tmp/MNI2Conform_nonl_fsl_field --absout --out=$M2M_DIR/toMNI/Conform2MNI_nonl
# note: after inversion, MNI2Conform_nonl_fsl_field defacto stores the warp field from
# conformed space to "MNI FSL mm space", so that concatenating it with the MNI FSL mm space to 
# MNI mm world space trafo (postmat.mat) gives a trafo from conformed mm world space to MNI mm world space


#############################
# create masks and files for visual inspection 
#############################
e echo 'writing files for visual control...'

# gm and wm masks in MNI space (non-binarized)
e subject2mni -i $M2M_DIR/wm_fromMesh.nii.gz -m $M2M_DIR -o $M2M_DIR/toMNI/wm # ending "MNI" will be automatically added
e subject2mni -i $M2M_DIR/gm_fromMesh.nii.gz -m $M2M_DIR -o $M2M_DIR/toMNI/gm

e subject2mni -i $M2M_DIR/T1fs_nu_conform.nii.gz -m $M2M_DIR -o $M2M_DIR/toMNI/T1fs_nu_nonlin
e subject2mni -i $M2M_DIR/T1fs_nu_conform.nii.gz -m $M2M_DIR -o $M2M_DIR/toMNI/T1fs_nu_12DOF -t 12dof 


#############################
# converting electrode csv files from MNI to subject space
#############################
e echo 'converting electrode csv files...'
e mni2subject_coords -m $M2M_DIR -s $SIMNIBSDIR/resources/ElectrodeCaps_MNI/EEG10-10_UI_Jurak_2007.csv -o $EEG_DIR/EEG10-10_UI_Jurak_2007.csv
