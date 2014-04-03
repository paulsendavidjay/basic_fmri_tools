#!/bin/sh

if [ $# -lt 4 ]; then
    echo "MUST SPECIFY SUBJECTID, VISIT, NUM_RUNS, \"RUN_LIST\""
    exit 0
fi


# DEFINE SUBJECT INFORMATION
SUBJECT_ID="$1"
VISIT="$2"
NUM_RUNS="$3"
RUN_LIST="$4"

RUN_LIST="${RUN_LIST%\"}" # remove end quote
RUN_LIST="${RUN_LIST#\"}" # remove front quote

echo $RUN_LIST
# SET LINE INFO TO $RUN
set -- $RUN_LIST
RUN1_NUM=$1
RUN2_NUM=$2
RUN3_NUM=$3
RUN4_NUM=$4

# SET UP DIRECTORY AND TEMPLATE FILE VARIABLES

if [ ${VISIT} -eq 1 ]; then
	DATA_DIR=data
	RUN_DIR=/Volumes/Phillips/bars/${DATA_DIR}/${SUBJECT_ID} # FOR PLACING FSF
else
	DATA_DIR=data_t${VISIT}
	RUN_DIR=/Volumes/Phillips/bars/${DATA_DIR}/${SUBJECT_ID} # FOR PLACING FSF
fi

# determine standard brain to use (depends on system)
if [ -e "/opt/ni_tools/fsl/data/standard/MNI152_T1_2mm_brain.nii.gz" ]; then
	STANDARD="/opt/ni_tools/fsl/data/standard/MNI152_T1_2mm_brain.nii.gz"
else
	STANDARD="/usr/local/fsl/data/standard/MNI152_T1_2mm_brain.nii.gz"
fi


TEMPLATEDIR=/Volumes/Phillips/bars/code_t2/FEAT_TEMPLATES

LVL1_MODEL_DIR=FEAT_MODEL2_whole_trial
LVL2_MODEL_DIR=FEAT_MODEL2_LVL2_whole_trial
FSF_NAME=${RUN_DIR}/FEAT_MODEL2_whole_trial.fsf # WHERE THE FSF COPY WILL BE MADE

if [ $NUM_RUNS -eq 2 ]; then


	LVL2_TEMPLATE=FEAT_LVL2_2_RUNS_whole_trial.fsf # WHAT WILL BE COPIED
	cd ${TEMPLATEDIR}
	for i in ${LVL2_TEMPLATE}; do
	 sed -e 's@DATA_DIR@'$DATA_DIR'@g' \
		-e 's@_STANDARD_@'$STANDARD'@g' \
		-e 's@LVL1_MODEL_DIR@'$LVL1_MODEL_DIR'@g' \
		-e 's@LVL2_MODEL_DIR@'$LVL2_MODEL_DIR'@g' \
		-e 's@SUBJECT_ID@'$SUBJECT_ID'@g' \
		-e 's@RUN1_NUM@'$RUN1_NUM'@g' \
		-e 's@RUN2_NUM@'$RUN2_NUM'@g'  <$i> ${FSF_NAME}
	done

	if [ -e ${RUN_DIR}/${LVL2_MODEL_DIR}.gfeat ]; then
		rm -R ${RUN_DIR}/${LVL2_MODEL_DIR}.gfeat
	fi

	# RUN FEAT SCRIPT
	feat ${FSF_NAME}


elif [ $NUM_RUNS -eq 3 ]; then


	LVL2_TEMPLATE=FEAT_LVL2_3_RUNS_whole_trial.fsf # WHAT WILL BE COPIED
	cd ${TEMPLATEDIR}
	for i in ${LVL2_TEMPLATE}; do
	 sed -e 's@DATA_DIR@'$DATA_DIR'@g' \
		-e 's@_STANDARD_@'$STANDARD'@g' \
		-e 's@LVL1_MODEL_DIR@'$LVL1_MODEL_DIR'@g' \
		-e 's@LVL2_MODEL_DIR@'$LVL2_MODEL_DIR'@g' \
		-e 's@SUBJECT_ID@'$SUBJECT_ID'@g' \
		-e 's@RUN1_NUM@'$RUN1_NUM'@g' \
		-e 's@RUN2_NUM@'$RUN2_NUM'@g' \
		-e 's@RUN3_NUM@'$RUN3_NUM'@g'  <$i> ${FSF_NAME}
	done
	
	if [ -e ${RUN_DIR}/${LVL2_MODEL_DIR}.gfeat ]; then
		rm -R ${RUN_DIR}/${LVL2_MODEL_DIR}.gfeat
	fi
	# RUN FEAT SCRIPT
	feat ${FSF_NAME}


elif [ $NUM_RUNS -eq 4 ]; then
	LVL2_TEMPLATE=FEAT_LVL2_4_RUNS_whole_trial.fsf # WHAT WILL BE COPIED
	cd ${TEMPLATEDIR}
	for i in ${LVL2_TEMPLATE}; do
	 sed -e 's@DATA_DIR@'$DATA_DIR'@g' \
		-e 's@_STANDARD_@'$STANDARD'@g' \
		-e 's@LVL1_MODEL_DIR@'$LVL1_MODEL_DIR'@g' \
		-e 's@LVL2_MODEL_DIR@'$LVL2_MODEL_DIR'@g' \
		-e 's@SUBJECT_ID@'$SUBJECT_ID'@g' <$i> ${FSF_NAME} # DIRECTORIES FOR SUBJECTS WITH FOUR RUNS HARD CODED
	done

	if [ -e ${RUN_DIR}/${LVL2_MODEL_DIR}.gfeat ]; then
		rm -R ${RUN_DIR}/${LVL2_MODEL_DIR}.gfeat
	fi

	# RUN FEAT SCRIPT
	feat ${FSF_NAME}
fi



