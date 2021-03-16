##################################################
# functions called by the other parts of dwi2cond
##################################################


e() {
  # only stderr is stored to logfile
  # does not work for all commands, due to redirection (e.g. do not use for cat) and the fact that
  # the command will be executed in a subshell (e.g. do not use for cd)
  echo '' >> $LOGFILE
  echo "$@" >> $LOGFILE
  "$@" 3>&1 1>&2 2>&3 3>&- | tee -a $LOGFILE

  ERRCODE=${PIPESTATUS[0]}
  if [ $ERRCODE -gt 0 ] && $EXIT_ON_ERROR; then
	echo 'exiting dwi2cond after error'
	echo 'error code' $ERRCODE >> $LOGFILE
	echo 'exiting dwi2cond' >> $LOGFILE
	echo '</pre></BODY></HTML>' >> $LOGFILE
	exit 1
  fi
}


e2() {
  # logs error code rather than error text to file; also works for cd and cat
  echo '' >> $LOGFILE
  echo "$@" >> $LOGFILE
  "$@"

  ERRCODE=$?
  if [ $ERRCODE -gt 0 ] && $EXIT_ON_ERROR; then
	echo 'exiting dwi2cond after error'
	echo 'error code' $ERRCODE >> $LOGFILE
	echo 'exiting dwi2cond' >> $LOGFILE
	echo '</pre></BODY></HTML>' >> $LOGFILE
	exit 1
  fi
}


APPLY_EDDY() {
# apply fsl eddy
#
# the following files will be created:
# DWI_corr: eddy-corrected file
# DWI_mean: mean image of DWI_corr
# DWI_cv: coefficient-of-variation image of DWI_corr
# 
# USAGE:
# without fsl topup results:
# APPLY_EDDY DWIraw DWIbvecs DWIbvals DWI_acqp DWI_index nodif_brain_mask
#
# with fsl topup results:
# APPLY_EDDY DWIraw DWIbvecs DWIbvals DWI_acqp DWI_index nodif_brain_mask_UNDIST topup_res

  local DWIRAW DWI_BVECS DWI_BVALS DWI_ACQP DWI_INDEX NODIF_BRAINMASK TOPUP_RES OPTSTR

  DWIRAW=$1
  DWI_BVECS=$2
  DWI_BVALS=$3
  DWI_ACQP=$4
  DWI_INDEX=$5
  NODIF_BRAINMASK=$6
  TOPUP_RES=$7

  e2 echo 'running fsl eddy ...'
  e2 echo "using raw dwi data" $DWIRAW
  e2 echo "using b-vectors" $DWI_BVECS
  e2 echo "using b-values" $DWI_BVALS
  e2 echo "using acquisition parameters" $DWI_ACQP
  e2 echo "using dwi index file" $DWI_INDEX
  e2 echo "using brain mask" $NODIF_BRAINMASK

  OPTSTR=''
  if [ -n "$TOPUP_RES" ]; then
     e2 echo ' using fsl topup results...' $TOPUP_RES
     OPTSTR='--topup='$TOPUP_RES
  fi

  e $EDDYBINARY --imain=$DWIRAW --mask=$NODIF_BRAINMASK --bvecs=$DWI_BVECS --bvals=$DWI_BVALS \
       --out=DWI_corr --acqp=$DWI_ACQP --index=$DWI_INDEX $OPTSTR $EDDYOPT --verbose
}


APPLY_FMCORR() {
# prepare FM for distortion correction
#
# result files (needed by EDDY_CORRECT):
# DWI_warp
# nodif_brain_mask_UNDIST
# nodif_brain_UNDIST
#
# USAGE: APPLY_FMCORR nodif_brain nodif_brainmask FMMAGraw FMPHASEraw fieldmap_TEdiff DWI_dwelltime warp-direction dodenoise

  local NODIF_BRAIN NODIF_BMASK FMMAGRAW FMPHASERAW MFMTEDIFF DWIDWELL UDIR DENOISE v OLD_PWD

  NODIF_BRAIN=$1
  NODIF_BMASK=$2
  FMMAGRAW=$3
  FMPHASERAW=$4
  FMTEDIFF=$5
  DWIDWELL=$6
  UDIR=$7
  DENOISE=$8
 
  e2 echo "preparing fieldmap to correct for B0 inhomogeneity..."
  e2 echo "using b=0 (brain)" $NODIF_BRAIN
  e2 echo "using brain mask" $NODIF_BMASK
  e2 echo "using field map magnitude" $FMMAGRAW
  e2 echo "using field map phase" $FMPHASERAW
  e2 echo "Fieldmap settings:"
  e2 echo "fieldmap_TEdiff " $FMTEDIFF
  e2 echo "DWI_DWELL" $DWIDWELL
  e2 echo "warping direction" $UDIR
  e2 echo "median filtering for denoising" $DENOISE

  # extract first subvolume of magnitude image
  e fslroi FMMAGraw FM_MAG 0 1 
  e bet FM_MAG FM_MAG_brain -m

  if [ ! "$FMTEDIFF" = -1 ]; then
    e2 echo 'Scaling fieldmap to [rad/s] (assuming Siemens gre field-map!)'
    # rescale the delta phase to be between 0 and 2*pi; starts out at -4096..+4096
    # (Siemens specific)
    e fsl_prepare_fieldmap SIEMENS $FMPHASERAW FM_MAG_brain FM_PHS_RADS $FMTEDIFF
  else
    e2 echo 'Assuming that fieldmap is already scaled in [rad/s]'
    e fslmaths $FMPHASERAW FM_PHS_RADS -odt float
  fi;

  # apply some de-noising to the phase map
  if $DENOISE; then
    e2 echo 'denoising phase map'
    e fslmaths FM_PHS_RADS -fmedian FM_PHS_RADS
  fi

  # get rid of constant signal offsets
  v=`fslstats FM_PHS_RADS -k FM_MAG_brain_mask -P 50`
  e fslmaths FM_PHS_RADS -sub $v FM_PHS_RADS 

  # creating distorted FM magnitude image as registration target
  # remark: not using sigloss, as DWI is not a gradient echo EPI
  DWIDWELL=`echo "scale=5; $DWIDWELL/1000" | bc`
  e fugue -i FM_MAG_brain --loadfmap=FM_PHS_RADS --mask=FM_MAG_brain_mask --dwell=$DWIDWELL \
  -w FM_MAG_brain_DIST --nokspace --unwarpdir=$UDIR

  # register DWI b=0 image to distorted FM magnitude image
  e flirt -in $NODIF_BRAIN -ref FM_MAG_brain_DIST -omat nodif2FM.mat -o nodif_brain_FMspace \
  -dof 6 -cost mutualinfo -searchcost mutualinfo
  e convert_xfm -omat FM2nodif.mat -inverse nodif2FM.mat
  
  e flirt -in FM_PHS_RADS -ref $NODIF_BRAIN -init FM2nodif.mat -applyxfm -out FM_PHS_RADS_DWIspace
  e flirt -in FM_MAG -ref $NODIF_BRAIN -init FM2nodif.mat -applyxfm -out FM_MAG_DWIspace
  e flirt -in FM_MAG_brain_mask -ref $NODIF_BRAIN -init FM2nodif.mat -applyxfm \
  -out FM_MAG_brain_mask_DWIspace 
  e fslmaths FM_MAG_brain_mask_DWIspace -thr 0.5 -bin FM_MAG_brain_mask_DWIspace -odt float
  
  # apply field map to DWI b=0 image
  e2 echo 'getting voxel shift map for DWI data'
  e fugue --loadfmap=FM_PHS_RADS_DWIspace --dwell=$DWIDWELL -i $NODIF_BRAIN -u nodif_brain_UNDIST \
  --unwarpdir=$UDIR --saveshift=DWI_shift --mask=FM_MAG_brain_mask_DWIspace
  e convertwarp -s DWI_shift -o DWI_warp -r $NODIF_BRAIN --shiftdir=$UDIR

  e applywarp -i $NODIF_BMASK -o nodif_brain_mask_UNDIST -w DWI_warp -r $NODIF_BRAIN --abs --interp=sinc
  e fslmaths nodif_brain_mask_UNDIST -mul FM_MAG_brain_mask_DWIspace -bin nodif_brain_mask_UNDIST
  e fslmaths nodif_brain_UNDIST -mul nodif_brain_mask_UNDIST nodif_brain_UNDIST
  RETURNVAL=$?

  rm FM_MAG* FM_PHS_* nodif_brain_FMspace.*

  return $RETURNVAL
}


APPLY_TOPUP() {
# apply fsl topup
# 
# result files (used by APPLY_EDDY):
# nodif_brain_mask_UNDIST
# topup_res_fieldcoef
# topup_res_movpar
#
# USAGE: APPLY_TOPUP nodif b0_reversedphase DWI_acqp config-file.cnf

  nodif=$1
  DWIrevphase=$2
  DWIacqp=$3
  conffile=$4

  e2 echo 'running fsl topup ...'
  e2 echo "using b=0 file" ${nodif}
  e2 echo "using b=0 file with reversed phase" ${DWIrevphase}
  e2 echo "using topup acqp file" ${DWIacqp}
  e2 echo "using topup config file" ${conffile}
  
  # concatenate nodif and nodif with reversed phase, then run topup
  e fslmerge -t nodif_topup ${nodif} ${DWIrevphase}
  e topup --imain=nodif_topup --datain=${DWIacqp} --config=${conffile} --out=topup_res \
  --iout=topup_iout --logout=topup.settings --verbose

  # get undistorted brain mask for fsl eddy
  e fslroi topup_iout nodif_UNDIST 0 1
  e bet nodif_UNDIST tmpBr -f 0.2 -m
  e fslmaths tmpBr nodif_brain_UNDIST
  e fslmaths tmpBr_mask nodif_brain_mask_UNDIST
  RETURNVAL=$?

  # tidy up
  e rm tmpBr* topup_iout.* nodif_topup.* nodif_UNDIST.*
  return $RETURNVAL
}


EDDY_CORRECT () {
# replacement of FSL eddy_correct
#  
# the following file will be created:
# DWI_corr: eddy-corrected file
# 
# USAGE:
# without fieldmap results
# EDDY_CORRECT DWIraw DWIbvals nodif nodif_brain_mask
#
# with fieldmap results
# EDDY_CORRECT DWIraw DWIbvals nodif nodif_brain_mask warpfile

  local DWIRAW DWI_BVALS NODIF NODIF_BRAINMASK WARPFILE FMCORR REGMAT REGMAT_ARR i BVALUES COUNTER NUMSTR FNAME_TMP FNAME_DWIREG
  
  DWIRAW=$1
  DWI_BVALS=$2
  NODIF=$3
  NODIF_BRAINMASK=$4

  FMCORR=false
  if [ $# -gt 5 ]; then
    WARPFILE=$5
    e2 echo "B0 unwarping:" ${WARPFILE}
    FMCORR=true
  fi

  e2 echo "eddy-current correction: pass 1"
  e GET_MEAN $DWIRAW DWIraw_mean $DWI_BVALS
  e mcflirt -in $DWIRAW -o DWI_pass1_corr -dof 6 -reffile DWIraw_mean
  e GET_MEAN DWI_pass1_corr DWI_pass1_mean $DWI_BVALS

  e2 echo "eddy-current correction: pass 2"
  e mcflirt -in $DWIRAW -o DWI_corr -dof 12 -sinc_final -reffile DWI_pass1_mean -mats
  e GET_MEAN DWI_corr DWI_corr_mean $DWI_BVALS # mean image for registration to nodif

  # get registration from mean DWI image to nodif, and concatenate with motion parameters
  e flirt -in DWI_corr_mean -ref $NODIF -nosearch -cost mutualinfo -interp sinc -omat meanDWI2nodif.mat -o meanDWI2nodif_tst
  REGMAT=$(ls DWI_corr.mat/MAT_*)
  REGMAT_ARR=( $REGMAT )
  for i in "${REGMAT_ARR[@]}" ; do
    e convert_xfm -omat ${i} -concat meanDWI2nodif.mat ${i}
  done
  e rm meanDWI2nodif_tst.*
  e rm meanDWI2nodif.mat

  # get registrations for b=0 images, and replace corresponding transformation files in DWI_corr.mat/
  e2 echo "registering b=0 images to nodif"
  e mcflirt -in $DWIRAW -o DWI_b0 -dof 6 -reffile $NODIF -mats
  
  BVALUES=$(<$DWI_BVALS)
  COUNTER=0
  for i in $BVALUES; do
     if [ $i -eq 0 ]; then
	NUMSTR=`printf "%04d\n" $COUNTER`
        e cp DWI_b0.mat/MAT_$NUMSTR DWI_corr.mat/MAT_$NUMSTR
     fi
     let COUNTER++
  done
  e rm DWI_b0.mat/*
  e rmdir DWI_b0.mat
  e rm DWI_b0.*

  # apply registration to DWIraw to create final DWI_corr
  e fslsplit $DWIRAW DWI_tmp
  FNAME_TMP=(`imglob DWI_tmp????.*`)

  if ! $FMCORR; then
    e2 echo 'applying registration to' ${FNAME_TMP[0]} '...'
    COUNTER=0;
    for i in "${FNAME_TMP[@]}" ; do
      echo $i ${REGMAT_ARR[$COUNTER]}
      e flirt -in $i -ref $NODIF -applyxfm -paddingsize 1 -interp sinc -init ${REGMAT_ARR[$COUNTER]} -o ${i}_reg
      let COUNTER++
    done

  else
    e2 echo 'applying warping+registration to' ${FNAME_TMP[0]} '...'
    COUNTER=0;
    for i in "${FNAME_TMP[@]}" ; do
      echo $i ${REGMAT_ARR[$COUNTER]}
      e applywarp -i $i -o ${i}_reg --premat=${REGMAT_ARR[$COUNTER]} -w ${WARPFILE} -r $NODIF --abs --interp=sinc
      let COUNTER++
    done
  fi

  FNAME_DWIREG=(`imglob DWI_tmp????_reg.*`)
  e fslmerge -t DWI_corr ${FNAME_DWIREG[@]} # final corrected image (unwarped + coregistered to nodif)
  RETURNVAL=$?

  # final mean image (for control only)
  FNAME_DWONLY=()
  COUNTER=0
  for i in $BVALUES; do
     if [ $i -gt 0 ]; then
	FNAME_DWONLY+=(${FNAME_DWIREG[$COUNTER]})
     fi
     let COUNTER++
  done
  e fslmerge -t DWI_corr_mean ${FNAME_DWONLY[@]}
  e fslmaths DWI_corr_mean -Tmean DWI_corr_mean

  rm DWI_pass1_*
  rm DWI_tmp????.*
  rm DWI_tmp????_reg.*
  rm DWI_corr.mat/*
  rmdir DWI_corr.mat

  return $RETURNVAL
}


GET_MEAN () {
# get mean image of B>0 images
# USAGE: GET_MEAN DWIdata fname_out bvalfile
  local INPUTFILE OUTPUTFILE BVALFILE TMPDATA FNAME_TMP BVALUES FNAME_DWONLY COUNTER i
  INPUTFILE=$1
  OUTPUTFILE=$2
  BVALFILE=$3

  TMPDATA=`remove_ext ${INPUTFILE}`
  e2 echo "getting mean of" $TMPDATA 
  e2 echo "reading b-vals from" ${BVALFILE}

  e fslsplit $TMPDATA DWI_tmp
  FNAME_TMP=(`imglob DWI_tmp????.*`)
  BVALUES=$(<${BVALFILE})
  FNAME_DWONLY=()
  COUNTER=0
  for i in $BVALUES; do
     if [ $i -gt 0 ]; then
	FNAME_DWONLY+=(${FNAME_TMP[$COUNTER]})
     fi
     let COUNTER++
  done
  e fslmerge -t $OUTPUTFILE ${FNAME_DWONLY[@]}
  e fslmaths $OUTPUTFILE -Tmean $OUTPUTFILE

  e2 echo "written file:" $OUTPUTFILE
  rm DWI_tmp????.*

  return $RETURNVAL
}



SAY() { 
  echo " "  | tee -a $LOGFILE
  echo "==============================" | tee -a $LOGFILE
  echo "==> $1 ..." | tee -a "$LOGFILE"
  echo "==============================" | tee -a $LOGFILE
}


WRITE_DWIACQP() {
# create files with acquistion parameters for fsl topup and eddy
#
# the following files will be created:
# ${basenameOut}_index
# ${basenameOut}_acqp
# see FSL docu for further information
#
# USAGE: WRITE_DWIACQP DWIrawdata basenameOut phase-direction readouttime
  local BASENAME PHASEDIR READOUTTIME
  DWIRAW=$1
  BASENAME=$2
  PHASEDIR=$3
  READOUTTIME=$4

  e2 echo 'writing acquistion parameters for fsl topup and/or eddy...'
  e2 echo 'Phase direction :' $PHASEDIR
  e2 echo 'Readout time: '$READOUTTIME' [s]'

  case "$PHASEDIR" in
      "x") PHASEVEC=' 1 0 0 '; PHASEVECREV='-1 0 0 '; NPHASE=`fslval $DWIRAW dim1| tr -d ' '`;;
     "-x") PHASEVEC='-1 0 0 '; PHASEVECREV=' 1 0 0 '; NPHASE=`fslval $DWIRAW dim1| tr -d ' '`;;
      "y") PHASEVEC='0  1 0 '; PHASEVECREV='0 -1 0 '; NPHASE=`fslval $DWIRAW dim2| tr -d ' '`;;
     "-y") PHASEVEC='0 -1 0 '; PHASEVECREV='0  1 0 '; NPHASE=`fslval $DWIRAW dim2| tr -d ' '`;;
      "z") PHASEVEC='0 0  1 '; PHASEVECREV='0 0 -1 '; NPHASE=`fslval $DWIRAW dim3| tr -d ' '`;;
     "-z") PHASEVEC='0 0 -1 '; PHASEVECREV='0 0  1 '; NPHASE=`fslval $DWIRAW dim3| tr -d ' '`;;
        *) echo "ERROR: Phase encoding direction "$PHASEDIR" unclear."; exit 1;;
  esac

  rm -f ${BASENAME}_acqp
  echo $PHASEVEC $READOUTTIME>>${BASENAME}_acqp
  echo $PHASEVECREV $READOUTTIME>>${BASENAME}_acqp

  NVOL=`fslval $DWIRAW dim4`
  rm -f ${BASENAME}_index
  for ((i=1; i<=NVOL; i++)); do
    printf '1 '>>${BASENAME}_index
  done
}

