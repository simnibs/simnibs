#!/bin/bash

# this script maps a few coil and electrode positions from MNI to subject space,
# prepares a simstruct and runs simnibs
# 
# INPUT: it needs the subject ID as input, and has to be called from within
# the directory containing the head mesh, the m2m_subid folder etc
#
# A. Thielscher, 2018

if [ $# -lt 1 ]; then
        echo "ERROR: this script needs the subject ID as input"; exit 1;
fi

SUBID=$1

pythoninstance=$SIMNIBSDIR/miniconda2/envs/simnibs_env/bin/python

TESTFILEDIR=$SIMNIBSDIR/Python_modules/src/simnibs/tests/simulation

PREPSCRIPT=$TESTFILEDIR/test_simnibsrun_prep_simstruct.py
PROTO_POS_CSV=$TESTFILEDIR/test_simnibsrun_positions.csv
PROTO_SIMSTRUCT=$TESTFILEDIR/test_simnibsrun_simstruct.mat

# map positions from MNI to subject space
mni2subject_coords -m m2m_${SUBID} -s $PROTO_POS_CSV -o testing_pos_${SUBID}.csv

# prepare simstruct (writes ${SUBID}_test_simstruct.mat)
$pythoninstance $PREPSCRIPT $SUBID testing_pos_${SUBID}.csv $PROTO_SIMSTRUCT

# run simulations
simnibs ${SUBID}_test_simstruct.mat

