#!/bin/bash

export SUBJECTS_DIR=/mnt/NAS_Progetti/PSP-MSA/daPellecchia/Bids_Salerno/derivatives/freesurfer
#subjects_list=$(find /mnt/NAS_Progetti/PSP-MSA/daPellecchia/Bids_Salerno/derivatives/freesurfer -maxdepth 1 -type d -name 'sub-HC*' | sort)
subject_list=$(find /mnt/NAS_Progetti/PSP-MSA/daPellecchia/Bids_Salerno/derivatives/freesurfer -maxdepth 1 -type d -name 'sub-*_ses-T0' | sort)
nthreads=16

for structure in brainstem; # {thalamus,hippo-amygdala} 
do
    for folder in ${subject_list};
    do 
        #if [[ ! "$folder" == *"long"* ]];then # exluding the longitudinal subjects, get only the base_
            #echo $folder
            #
            #echo $subj_id
            #segment_subregions ${structure} --long-base ${subj_id} --debug --threads ${nthreads}  # IF LONGITUDINAL DATA USE --long-base
            #segment_subregions ${structure} --cross ${subj_id} --debug --threads ${nthreads} # IF cross I ran this for the controls
        if [[ ! "$folder" == *"sub-HC"* ]];then  # for the cross study get the patients (the HC have been already done) 
            subj_id=$(basename $folder)
            echo $subj_id
            segment_subregions ${structure} --cross ${subj_id} --debug --threads ${nthreads} # run cross forthe patient sonly cross sectional (t0 all subejcts indluing also the one without the FUP)
        fi 

    done 
done 
