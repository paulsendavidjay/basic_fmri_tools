#!/bin/sh

SUBJECT_ID=$1
VISIT=$2
RUN_NUM=$3

if [ $# -lt 3 ]; then
    echo "MUST SPECIFY SUBJECTID, VISIT, RUN_NUM"
    exit 0
fi


# SET UP DIRECTORY AND TEMPLATE FILE VARIABLES
if [ ${VISIT} -eq 1 ]; then
	DATA_DIR=data
	RUN_DIR=/Volumes/Phillips/bars/${DATA_DIR}/${SUBJECT_ID}/bars_run${RUN_NUM}
	if [ ! -e ${RUN_DIR}/nfswkmtd_${SUBJECT_ID}_5.nii.gz ]; then
    	if [ ! -e ${RUN_DIR}/nfswkmtd_functional.nii.gz ]; then
			echo "${SUBJECT_ID} data not found" >> FEAT_errs_nodata.txt
			echo "${SUBJECT_ID} data not found"
			exit 0
		else
			FUNC_DATA=nfswkmtd_functional.nii.gz
		fi
	else
		FUNC_DATA=nfswkmtd_${SUBJECT_ID}_5.nii.gz
	fi
else
	DATA_DIR=data_t${VISIT}
	RUN_DIR=/Volumes/Phillips/bars/${DATA_DIR}/${SUBJECT_ID}/bars_run${RUN_NUM}
	if [ ! -e ${RUN_DIR}/nfswkmtd_functional.nii.gz ]; then
    	echo "${SUBJECT_ID} data not found" >> FEAT_errs_nodata.txt
	    echo "${SUBJECT_ID} data not found"
	    exit 0
	else
		FUNC_DATA=nfswkmtd_functional.nii.gz
	fi
fi



# skip analyis or remove old feat directory
if [ -e ${RUN_DIR}/FEAT_MODEL2_whole_trial.feat/stats/cope5.nii.gz ]; then
	echo "Already ran"
	exit
else
	rm -R ${RUN_DIR}/FEAT_MODEL2_whole_trial.feat
fi

# determine standard brain to use (depends on system)
if [ -e "/opt/ni_tools/fsl/data/standard/MNI152_T1_2mm_brain.nii.gz" ]; then
	STANDARD="/opt/ni_tools/fsl/data/standard/MNI152_T1_2mm_brain.nii.gz"
else
	STANDARD="/usr/local/fsl/data/standard/MNI152_T1_2mm_brain.nii.gz"
fi


TEMPLATEDIR=/Volumes/Phillips/bars/code_t2/FEAT_TEMPLATES
FEAT_OUTPUT_DIR=${RUN_DIR}/FEAT_MODEL2_whole_trial.feat # WHERE THE FEAT OUTPUT WILL GO

# TEMPLATE AND SUBJECT/RUN FSF 
LVL1_TEMPLATE=FEAT_LVL1_SUBJECT_ID_RUN_NUM_v5.98_whole_trial.fsf # WHAT WILL BE COPIED
FSF_NAME=${RUN_DIR}/FEAT_MODEL2_whole_trial.fsf # WHERE THE FSF COPY WILL BE MADE




cd ${TEMPLATEDIR}
for i in ${LVL1_TEMPLATE}; do
 sed -e 's@SUBJECT_ID@'$SUBJECT_ID'@g' \
 -e 's@RUN_NUM@'$RUN_NUM'@g' \
 -e 's@FUNC_DATA@'$FUNC_DATA'@g' \
 -e 's@_STANDARD_@'$STANDARD'@g' \
 -e 's@DATA_DIR@'$DATA_DIR'@g' <$i> ${FSF_NAME}
done


# RUN FEAT SCRIPT
feat ${FSF_NAME}

# CHECK TO MAKE SURE FEAT COMPLETED, REMOVE LARGE UNNECESSARY FILES
if [ -e ${FEAT_OUTPUT_DIR}/stats/corrections.nii.gz ]; then
	rm ${FEAT_OUTPUT_DIR}/stats/corrections.nii.gz # approx 120Mb
	rm ${FEAT_OUTPUT_DIR}/stats/res4d.nii.gz # approx 90Mb
else
    echo "${SUBJECT_ID} ${VISIT} ${RUN_NUM} FEAT DID NOT COMPLETE" >> FEAT_errs_incomplete.txt
    echo "${SUBJECT_ID} ${VISIT} ${RUN_NUM} FEAT DID NOT COMPLETE"
    exit 0
fi

# UPSAMPLE COPE AND VARCOPE IMAGES
#featregapply ${FEAT_OUTPUT_DIR}




