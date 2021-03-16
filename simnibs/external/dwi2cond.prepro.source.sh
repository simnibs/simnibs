
OLD_PWD=`pwd`

if $DTI_EXISTS; then
# using existing diffusion tensor (output from fsl dtifit)
###########################################################
  e echo "Copying DTI tensor"

  if [ ! -f $DWIDATA ]; then 
     e echo "ERROR: Input file '$DWIDATA' doesn't exist."
     exit 1
  fi

  if [ ! -d $DTIRAW_DIR ]; then e mkdir $DTIRAW_DIR; fi 
  e fslmaths $DWIDATA $DTIRAW_DIR/DTI_tensor
  e fslreorient2std $DTIRAW_DIR/DTI_tensor $DTIRAW_DIR/DTI_tensor
  e fslmaths $DTIRAW_DIR/DTI_tensor -tensor_decomp $DTIRAW_DIR/DTI
  e rm $DTIRAW_DIR/DTI.*

else
# preprocessing of DWI raw data
###############################
  
# check for raw input files
  for i in {$DWIDATA,$BVALS,$BVECS,$FMMAG,$FMPHASE,$REVPHASEDATA}; do
      if [ ! -f $i ]; then 
        echo "ERROR: Input file '$i' doesn't exist."
        exit 1
      fi
  done


# copy raw files to subdirectories
  if [ ! -d $RAW_DIR ]; then e mkdir $RAW_DIR; fi 
  e fslmaths $DWIDATA $RAW_DIR/DWIraw # copying with fslmaths ensures that file type is changed to .nii.gz
  e fslreorient2std $RAW_DIR/DWIraw $RAW_DIR/DWIraw # this prevents L/R flips of tensors for neurological input data
  e cp $BVALS $RAW_DIR/DWIbvals

  case `uname` in
    Linux)
	  e sed -i 's/\r//' $RAW_DIR/DWIbvals # remove any CR characters
    ;;
    Darwin)
	  e sed -i ''  's/\r//' $RAW_DIR/DWIbvals # remove any CR characters
    ;;
  esac
  e cp $BVECS $RAW_DIR/DWIbvecs

  if $FMCORR; then
    if [ ! -d $FM_DIR ]; then e mkdir $FM_DIR; fi 
    e fslmaths $FMMAG $FM_DIR/FMMAGraw
    e fslmaths $FMPHASE $FM_DIR/FMPHASEraw
    e fslreorient2std $FM_DIR/FMMAGraw $FM_DIR/FMMAGraw
    e fslreorient2std $FM_DIR/FMPHASEraw $FM_DIR/FMPHASEraw
  fi

  if $DOTOPUP; then
    if [ ! -d $TOPUP_DIR ]; then e mkdir $TOPUP_DIR; fi
    e cp "$DWI2CONDPATH/b02b0_nosubsamp.cnf" $TOPUP_DIR
    e fslmaths $REVPHASEDATA $TOPUP_DIR/DWIrevphase
    e fslreorient2std $TOPUP_DIR/DWIrevphase $TOPUP_DIR/DWIrevphase
    e mcflirt -in $TOPUP_DIR/DWIrevphase -o $TOPUP_DIR/DWIrevphase -dof 6
    e fslmaths $TOPUP_DIR/DWIrevphase -Tmean $TOPUP_DIR/DWIrevphase
  fi


# getting mean B=0 image and apply brain extraction
  OLD_PWD=`pwd`
  e2 cd $RAW_DIR
  e echo " "
  e echo "Preprocessing DWI data..."
  
  e echo "getting mean b0 image"
  e fslsplit DWIraw DWI_tmp
  FNAME_TMP=(`imglob DWI_tmp????.*`)
  BVALUES=$(<DWIbvals)
  FNAME_B0=()
  N_B0=0
  COUNTER=0
  for i in $BVALUES; do
     if [ $i -eq 0 ]; then
	FNAME_B0+=(${FNAME_TMP[$COUNTER]})
	let N_B0++
     fi
     let COUNTER++
  done
  e echo "found" $N_B0 "b0 volumes"
  e fslmerge -t nodif_all ${FNAME_B0[@]}
  e mcflirt -in nodif_all -o nodif -dof 6
  e fslmaths nodif -Tmean nodif
  e bet nodif nodif_brain -f 0.2 -m

  e rm DWI_tmp????.*
  e rm nodif_all.*
  
  # get SSE of raw data for QA
  e dtifit --data=DWIraw --bvecs=DWIbvecs --bvals=DWIbvals --out=DWIraw_SSE --sse --mask=nodif_brain_mask $DTIFITOPT
  e immv DWIraw_SSE_sse DWIraw_SSE
  e immv DWIraw_SSE_FA DWIraw_FA
  e rm DWIraw_SSE_*.*


# distortion correction and eddy current correction  
  if [ ! -d $DTIRAW_DIR ]; then e mkdir $DTIRAW_DIR; fi
  e cp $RAW_DIR/DWIbvals $DTIRAW_DIR # needed for dtifit, and not changed during eddy correction

  if $NOMOCO; then
    # no eddy current correction
    # --------------------------
    e echo "skipping eddy current correction ..."
    e imcp $RAW_DIR/DWIraw $DTIRAW_DIR/DWIforfit
    e cp $RAW_DIR/DWIbvecs $DTIRAW_DIR
    e imcp $RAW_DIR/nodif_brain_mask $DTIRAW_DIR
    e imcp $RAW_DIR/nodif_brain $DTIRAW_DIR


  elif $DOEDDY; then
    # FSL eddy and topup
    # --------------------------
    if [ ! -d $EDDY_DIR ]; then e mkdir $EDDY_DIR; fi
    e2 cd $EDDY_DIR
    # write file with acquistion parameters
    e WRITE_DWIACQP $RAW_DIR/DWIraw DWI ${PHASEDIR} ${READOUTTIME}
    
    TOPUP_STR=''
    if $DOTOPUP; then
      # run fsl topup (created files: nodif_brain_mask_UNDIST nodif_brain_UNDIST topup_res_fieldcoef topup_res_movpar)     
      if [ ! -d $TOPUP_DIR ]; then e mkdir $TOPUP_DIR; fi
      e2 cd $TOPUP_DIR
      TOPUP_STR="$TOPUP_DIR/topup_res"
      e APPLY_TOPUP $RAW_DIR/nodif DWIrevphase $EDDY_DIR/DWI_acqp b02b0_nosubsamp.cnf  
      e imcp $TOPUP_DIR/nodif_brain_mask_UNDIST $DTIRAW_DIR/nodif_brain_mask
      e imcp $TOPUP_DIR/nodif_brain_UNDIST $DTIRAW_DIR/nodif_brain
      e2 cd $EDDY_DIR
    else
      e imcp $RAW_DIR/nodif_brain_mask $DTIRAW_DIR
      e imcp $RAW_DIR/nodif_brain $DTIRAW_DIR
    fi

    # run fsl eddy (created files: DWI_corr DWI_corr.eddy_rotated_bvecs ...)
    e APPLY_EDDY $RAW_DIR/DWIraw $RAW_DIR/DWIbvecs $RAW_DIR/DWIbvals DWI_acqp DWI_index $DTIRAW_DIR/nodif_brain_mask $TOPUP_STR
    e cp DWI_corr.eddy_rotated_bvecs $DTIRAW_DIR/DWIbvecs
    e imcp DWI_corr $DTIRAW_DIR/DWIforfit
    

  else
    # standard eddy current correction and field-map based distortion correction
    # --------------------------------------------------------------------------
    if [ ! -d $EDDYCOR_DIR ]; then e mkdir $EDDYCOR_DIR; fi
    e2 cd $EDDYCOR_DIR
    
    FM_STR=''
    if $FMCORR; then 
      # run fieldmap-based distortion correction (created files: nodif_brain_mask_UNDIST nodif_brain_UNDIST DWI_warp)
      e2 cd $FM_DIR
      e APPLY_FMCORR $RAW_DIR/nodif_brain $RAW_DIR/nodif_brain_mask FMMAGraw FMPHASEraw $FMTEDIFF $DWIDWELL $UDIR $DENOISE
      FM_STR="$FM_DIR/DWI_warp"
      e imcp $FM_DIR/nodif_brain_mask_UNDIST $DTIRAW_DIR/nodif_brain_mask
      e imcp $FM_DIR/nodif_brain_UNDIST $DTIRAW_DIR/nodif_brain
      e2 cd $EDDYCOR_DIR
    else
      e imcp $RAW_DIR/nodif_brain_mask $DTIRAW_DIR
      e imcp $RAW_DIR/nodif_brain $DTIRAW_DIR
    fi
    e cp $RAW_DIR/DWIbvecs $DTIRAW_DIR

    # run standard eddy current correction (created files: DWI_corr)
    e EDDY_CORRECT $RAW_DIR/DWIraw $RAW_DIR/DWIbvals $RAW_DIR/nodif $RAW_DIR/nodif_brain_mask $FM_STR
    e imcp DWI_corr $DTIRAW_DIR/DWIforfit
  fi


# fitting of diffusion tensors
  e echo "getting DTI tensors"
  e2 cd $DTIRAW_DIR
  e fslmaths DWIforfit -thr 0 DWIforfit
  e dtifit -k DWIforfit -m nodif_brain_mask -r DWIbvecs -b DWIbvals -o DTI --sse --save_tensor $DTIFITOPT
  e rm DTI_V?.* DTI_L?.* DTI_M?.* DTI_S0.*


fi # if $DTI_EXISTS; then ... else
e2 cd $OLD_PWD


