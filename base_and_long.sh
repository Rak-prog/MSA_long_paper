#!/bin/bash

num_threads=8

anat_folder=${1}
output_path="/home/riccardo/codici_progetti/Salerno/tmp_freesurfer"  
freesurfer_nas="/mnt/NAS_Progetti/PSP-MSA/daPellecchia/Bids_Salerno/derivatives/freesurfer"
SUBJECTS_DIR=${output_path}
export SUBJECTS_DIR

while read subj_path; do

        subj_id=$(basename $subj_path | cut -d "_" -f 1)
        session=$(basename $subj_path | cut -d "_" -f 2)
        directory=$(dirname $subj_path)
        echo $directory
        echo $subj_path
        echo $subj_id
        echo $session
        count=$(find "$subj_path" -maxdepth 1 -type d -name 'ses-*' | wc -l)

        echo $count

        #if [ $count == "2" ]; then
            # copy these file from the NAS ${subj_id}_ses-T0 ${subj_id}_ses-T1 # se non funziona prova con opzione L in rsync
            rsync -aHAX --progress ${freesurfer_nas}/${subj_id}_ses-T0 ${output_path}
            rsync -aHAX --progress ${freesurfer_nas}/${subj_id}_ses-T1 ${output_path}

            # run unbiased template
            recon-all -base ${subj_id}_ -tp ${subj_id}_ses-T0 -tp ${subj_id}_ses-T1 -all -nthreads ${num_threads} #-parallel

            # run long 
            recon-all -long ${subj_id}_ses-T0 ${subj_id}_ -all -nthreads ${num_threads} #-parallel
            recon-all -long ${subj_id}_ses-T1 ${subj_id}_ -all -nthreads ${num_threads} #-parallel

            # copy the base and the long results to NAS
            rsync -aHAX --progress ${output_path}/${subj_id}_ ${freesurfer_nas}/ > ./rsync_logfile/${subj_id}_unbiased_template_local2NAS.txt
            rsync -aHAX --progress ${output_path}/${subj_id}_ses-T0.long.${subj_id}_ ${freesurfer_nas}/ > ./rsync_logfile/${subj_id}_long_local2NAS_T0.txt
            rsync -aHAX --progress ${output_path}/${subj_id}_ses-T1.long.${subj_id}_ ${freesurfer_nas}/ > ./rsync_logfile/${subj_id}_long_local2NAS_T1.txt

            # remove the previosuly copied NAS subject and the base results from local to remove space 
            rm -rf ${output_path}/${subj_id}_
            rm -rf ${output_path}/${subj_id}_ses-T0.long.${subj_id}_
            rm -rf ${output_path}/${subj_id}_ses-T1.long.${subj_id}_
            rm -rf ${output_path}/${subj_id}_ses-T0
            rm -rf ${output_path}/${subj_id}_ses-T1
        #fi   

done < ${anat_folder}