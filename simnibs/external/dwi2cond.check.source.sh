
# show some resutls for quality control
############################################

  OLD_PWD=`pwd`
  e2 cd $D2C_DIR
  e cp $DWI2CONDPATH/d2c_check.opt .

  e2 cd $DTICONF_DIR
  FSLEYESBIN=$(which fsleyes)
  if [ ${#FSLEYESBIN} -gt 0 ]; then
    # check of final results
    e fsleyes T1_brain DTI_coregT1_FA -a 70 -dr 0.1 0.8 -cm red-yellow DTI_coregT1_V1 -a 0 -ot linevector &

    # check of sse of tensor fit and effect of distortion correction
    HLPSTR='DTI_FA_6dof_QA -dr 0.1 0.8 -cm red-yellow '
    if [ -f DTIraw_SSE_6dof_QA.nii.gz ]; then 
	HLPSTR=$HLPSTR'DTIraw_SSE_6dof_QA -dr 0.0 1.0 -cm hot '
    fi
    if [ -f DTI_SSE_6dof_QA.nii.gz ]; then 
	HLPSTR=$HLPSTR'DTI_SSE_6dof_QA -dr 0.0 1.0 -cm hot '
    fi
    e fsleyes T1_brain $HLPSTR T1_brainrim_QA -cm blue -a 40 &

  else
    # check of final results
    e fslview T1_brain DTI_coregT1_FA -t 0.7 -b 0.1,0.8 -l Red-Yellow DTI_coregT1_V1 -t 0 &

    # check of sse of tensor fit and effect of distortion correction
    HLPSTR='DTI_FA_6dof_QA -b 0.1,0.8 -l Red-Yellow '
    if [ -f DTIraw_SSE_6dof_QA.nii.gz ]; then 
	HLPSTR=$HLPSTR'DTIraw_SSE_6dof_QA -b 0.0,1.0 -l Hot '
    fi
    if [ -f DTI_SSE_6dof_QA.nii.gz ]; then 
	HLPSTR=$HLPSTR'DTI_SSE_6dof_QA -b 0.0,1.0 -l Hot '
    fi
    e fslview T1_brain $HLPSTR T1_brainrim_QA -l Blue -t 0.4 &
  fi

  e ${BINDIR}/gmsh $D2C_DIR/first_ev_for_check.msh $D2C_DIR/d2c_check.opt &

  e2 cd $OLD_PWD
