#!/bin/bash
export SUBJECTS_DIR=/mnt/NAS_Progetti/PSP-MSA/daPellecchia/Bids_Salerno/derivatives/freesurfer
group_qdec="/home/riccardo/codici_progetti/Salerno/qdec/two_timepoints_allMSA_long.qdec.table.dat"  # this is with the python script create_qdec_final.py 
group_qdec_cross_HC="/home/riccardo/codici_progetti/Salerno/qdec/qdec_controls.qdec.table.dat"
protocol_regressout=${1} # 'all', 'two', 'merge' all three protocols but un group is unbalnced only four subjects, two is without these four subjects and merged I merged the four to the group with the most similar feature the group with 22 MSA
do_fitting=0
do_plotting=1
correction_with_permutation=1
correction_with_MCZ=0

subj_id=${2} # if you want to visualize single subject result, option below visualzie set to 1

fsdg_files=(stats_results/MSA_between_group_MCI_${protocol_regressout}protocols_nosubj18nosubj45nosubj53/MSA_between_group_MCI.fsdg
stats_results/MSA_between_group_fenotipo_${protocol_regressout}protocols_nosubj18nosubj45nosubj53/MSA_between_group_fenotipo.fsdg 
stats_results/MSA_group_average_${protocol_regressout}protocols_nosubj18nosubj45nosubj53/MSA_average_group_age_sex.fsdg) 

#stats_results/MSAcerebellar_between_group_MCI_${protocol_regressout}protocols_nosubj18nosubj45nosubj53/MSAcerebellar_between_group_MCI.fsdg)

#fsdg_files=(stats_results/MSA_group_average_${protocol_regressout}protocols_nosubj18nosubj45nosubj53/MSA_average_group_age_protocol.fsdg)
#stats_results/MSAcerebellar_between_group_MCI_${protocol_regressout}protocols_nosubj18nosubj45nosubj53/MSAcerebellar_between_group_MCI.fsdg

output_stats="cluster" # voxel or cluster visualize the corrected ouput in freeview if are corrected at the cluster or at the voxel level

# control variables
detect_errors=0 # read all the recon-all log file and print thos folder with FINISHED WITH ERRORS
extract=0  # this is calling aparc2table amnd preparing the CT data from the long folder using the data in the qdec file (QDEC FILE MUST BE PROPERLY DEFINED with the python script create_qdec.py)
calc=0   # calcualte the percent change, rate of change of cortical thickness - TO DO use other measure like volumes, surface etcc..
visualize=0 # visualzie results of a selecte single subject. In this case is mandatory to pass an argument with the subject it 
dostat=1  # do concatenation of the variables of interests (i.e. symmetrize percent change) and run  a GLM t test only one gruoup -- TO IMPLEMENT betwenen group comparison
nthreads=16

# paramters for the permutation based correction are below:
NPERM=5000 #5000
TH=1.3 # this is 0.05 in -log10(p)
fwhm_value=10
long_parameter="spc"  # or rate or avg? check here https://surfer.nmr.mgh.harvard.edu/fswiki/LongitudinalTwoStageModel
low_p=3 # 3 is 0.001 and 1.3 is 0.05, 1.9 is 0.0125, 2 is 0.01
high_p=5
cluster_forming=3 # ,ust be 0.001 as in the paper of Greve and Fiscl 

low_p_corrected=1.3 # 3 is 0.001 and 1.3 is 0.05, 1.9 is 0.0125, 2 is 0.01
high_p_corrected=4

if [ $detect_errors == "1" ]; then
    files=$(find ${SUBJECTS_DIR}/sub-*_*/scripts/recon-all.log)
    echo "##### below the freesurfer cross, base or long that have crashed #####"
    #echo $files
    for f in $files;do  
        grep -i 'with error' $f; 
    done
fi 

if [ $extract == "1" ]; then

    for measure in {area,thicknessstd,thickness,volume,meancurv};
    do 
        for hemi in {lh,rh}; 
        #  --subjects sub-08_ses-T0.long.sub-08_ sub-09_ses-T0.long.sub-09_ sub-11_ses-T0.long.sub-11_
        do 
            # for cortical parcellation         
            aparcstats2table --qdec-long ${group_qdec} \
                            --meas ${measure} \
                            --hemi ${hemi} \
                            --parc aparc.DKTatlas \
                            --tablefile ./stats_results/${measure}_long_${hemi}.txt \
                            --parcid-only

            # cortical parcellation of controls only crosssections l(do not pass here qdec-long!!!)
            aparcstats2table --subjects sub-HC001_ses-T0 sub-HC002_ses-T0 sub-HC003_ses-T0 sub-HC004_ses-T0 sub-HC005_ses-T0 sub-HC006_ses-T0 sub-HC007_ses-T0 sub-HC008_ses-T0 sub-HC009_ses-T0 \
                                        sub-HC010_ses-T0 sub-HC011_ses-T0 sub-HC012_ses-T0 sub-HC013_ses-T0 sub-HC014_ses-T0 sub-HC015_ses-T0 sub-HC016_ses-T0 sub-HC017_ses-T0 sub-HC018_ses-T0 \
                                        sub-HC019_ses-T0 sub-HC020_ses-T0 sub-HC021_ses-T0 sub-HC022_ses-T0 sub-HC023_ses-T0 sub-HC024_ses-T0 sub-HC025_ses-T0 sub-HC026_ses-T0 sub-HC027_ses-T0 \
                                        sub-HC028_ses-T0 sub-HC029_ses-T0 sub-HC030_ses-T0 sub-HC031_ses-T0 sub-HC032_ses-T0 sub-HC033_ses-T0 sub-HC034_ses-T0 sub-HC035_ses-T0 sub-HC036_ses-T0 \
                                        sub-HC037_ses-T0 sub-HC038_ses-T0 sub-HC039_ses-T0 sub-HC040_ses-T0 sub-HC041_ses-T0 sub-HC042_ses-T0 sub-HC043_ses-T0 sub-HC044_ses-T0 sub-HC045_ses-T0 \
                                        sub-HC046_ses-T0 sub-HC047_ses-T0 sub-HC048_ses-T0 sub-HC049_ses-T0 sub-HC050_ses-T0 \
                            --meas ${measure} \
                            --hemi ${hemi} \
                            --parc aparc.DKTatlas \
                            --tablefile ./stats_results/${measure}_crossHC_${hemi}.txt \
                            --parcid-only

            # for subcortical aprcellation 
            if [[ "$measure" == *"volume"* ]]; then 
                asegstats2table --qdec-long ${group_qdec} \
                            --meas volume \
                            --all-segs \
                            --tablefile ./stats_results/${measure}_long_${hemi}_volumes.txt \
                            -v 
            fi 
        done
    done 
fi 

if [ $calc == "1" ]; then
    # this command calculate the rate of change (annualized) as I defined the years column in years for the --meas measure, each subject of the qdec table 
    # I am calculating per each subject (the other I left out for now at the moment these two seems the msot interesting for the longitudinal MSA study)
    # -do-rate compute the rate of change (thickening in mm/year) 
    # --do-spc compute the percent change with respect to the temporal average that is more stable and reliable than with respect to timpepoint 1)
    #long_mris_slopes --qdec ${group_qdec} --meas thickness --sd $SUBJECTS_DIR --hemi both --do-rate --do-avg --do-stack --do-spc --time years --qcache fsaverage --do-label  #--do-avg --do-rate --do-pc1 --do-stack --do-label
    long_stats_slopes --qdec ${group_qdec} --meas volume --stats=aseg.stats --sd $SUBJECTS_DIR --do-spc --time years --do-stack 
    long_stats_slopes --qdec ${group_qdec} --meas thickness --stats=lh.aparc.DKTatlas.stats --sd $SUBJECTS_DIR --do-spc --time years --do-stack 
    long_stats_slopes --qdec ${group_qdec} --meas thickness --stats=rh.aparc.DKTatlas.stats --sd $SUBJECTS_DIR --do-spc --time years --do-stack 
fi 

if [ $visualize == "1" ]; then 
    if [ -f $SUBJECTS_DIR/${subj_id}_/surf/lh.pial.T1 ]; then 
        cp $SUBJECTS_DIR/${subj_id}_/surf/lh.pial.T1 $SUBJECTS_DIR/${subj_id}_/surf/lh.pial
        cp $SUBJECTS_DIR/${subj_id}_/surf/rh.pial.T1 $SUBJECTS_DIR/${subj_id}_/surf/rh.pial
    fi 

    freeview -f $SUBJECTS_DIR/${subj_id}_/surf/lh.pial:overlay=$SUBJECTS_DIR/${subj_id}_/surf/lh.long.thickness-${long_parameter}.fwhm${fwhm_value}.mgh:overlay_threshold=1,3.5:overlay=$SUBJECTS_DIR/${subj_id}_/surf/lh.long.thickness-stack.mgh:annot=$SUBJECTS_DIR/${subj_id}_/label/lh.aparc.DKTatlas.annot:annot_outline=1 \
             -f $SUBJECTS_DIR/${subj_id}_/surf/rh.pial:overlay=$SUBJECTS_DIR/${subj_id}_/surf/rh.long.thickness-${long_parameter}.fwhm${fwhm_value}.mgh:overlay_threshold=1,3.5:overlay=$SUBJECTS_DIR/${subj_id}_/surf/rh.long.thickness-stack.mgh:annot=$SUBJECTS_DIR/${subj_id}_/label/rh.aparc.DKTatlas.annot:annot_outline=1 \
             --timecourse --colorscale --viewport 3d
            
    # visualize annualized CT for the subject passed as argument 
    #freeview -f $SUBJECTS_DIR/fsaverage/surf/lh.pial:overlay=$SUBJECTS_DIR/${subj_id}_/surf/lh.long.thickness-${long_parameter}.fwhm${fwhm_value}.fsaverage.mgh:overlay_threshold=2,5 \
    #        -f $SUBJECTS_DIR/fsaverage/surf/rh.pial:overlay=$SUBJECTS_DIR/${subj_id}_/surf/rh.long.thickness-${long_parameter}.fwhm${fwhm_value}.fsaverage.mgh:overlay_threshold=2,5 \
    #        --colorscale --viewport 3d
fi 


if [ $dostat == "1" ]; then 
    for fsdg_file in "${fsdg_files[@]}"; 
    do 
        # define here some important varibales to take for defining the GLM 
        if [[ "$fsdg_file" == *"MSAcerebellar_between_group_MCI"* ]]; then  # to do the MSA cerebellar only - THE OUTLIERS FOR MSA-c are removed already because they are all MSA-p
            if [ $protocol_regressout == "two" ]; then  # case cerebellar subejcts wihtout 1 that is in the protcolc P1 with only four subjects
                base_concat="/home/riccardo/codici_progetti/Salerno/stats_results/concatenate_file_path_MSAcerebellar_twoprotocols.txt" 
            else
                base_concat="/home/riccardo/codici_progetti/Salerno/stats_results/concatenate_file_path_MSAcerebellar.txt"
            fi 
        else   # otheriwse to the entire group of MSA withing group and between group comparison
            if [ $protocol_regressout == "two" ]; then # if I remove two subjects use this file instead 
                base_concat="/home/riccardo/codici_progetti/Salerno/stats_results/concatenate_file_path_noprotocol4subj_nooutliers.txt"  # file with the path prepared in python create_qdec_final.py 
            else 
                base_concat="/home/riccardo/codici_progetti/Salerno/stats_results/concatenate_file_path_nooutliers.txt" 
            fi 
        fi 

        basename ${fsdg_file} .fsdg
        if [[ "$fsdg_file" == *"average"* ]]; then 
            output_path="/home/riccardo/codici_progetti/Salerno/$(dirname ${fsdg_file})"
        elif [[ "$fsdg_file" == *"MCI"* ]]; then 
            output_path="/home/riccardo/codici_progetti/Salerno/$(dirname ${fsdg_file})"
        elif [[ "$fsdg_file" == *"fenotipo"* ]]; then 
            output_path="/home/riccardo/codici_progetti/Salerno/$(dirname ${fsdg_file})"
        fi 

        echo $output_path
        
        for hemi in {lh,rh}; do 

            name_concat=${hemi}.long.thickness-${long_parameter}.fwhm${fwhm_value}.fsaverage.mgh
           
            while IFS= read -r line; do
                echo "${line}/${name_concat}"
            done < ${base_concat} > ${output_path}/${hemi}_input_list.txt

            mri_concat --o ${output_path}/${hemi}.long.thickness-${long_parameter}.fwhm${fwhm_value}.fsaverage.4d.mgh \
                       --f ${output_path}/${hemi}_input_list.txt
            
            #mri_info ${output_path}/${hemi}.long.thickness-${long_parameter}.fwhm${fwhm_value}.fsaverage.4d.mgh 

            # before running the GLM check if the concatenation and the subject sin fsdg are in the correct order :) 
            cut -d'/' -f9 ${output_path}/${hemi}_input_list.txt | sort > ${output_path}/check1_temp.txt
            if [[ "$fsdg_file" == *"average"* ]];then
                cut -d' ' -f2 ${fsdg_file} | sort | tail -n +5 > ${output_path}/check2_temp.txt
            elif  [[ "$fsdg_file" == *"fenotipo"* ]] || [[ "$fsdg_file" == *"MCI"* ]]; then
                cut -d' ' -f2 ${fsdg_file} | sort | tail -n +6 > ${output_path}/check2_temp.txt
            fi 
            cut -d' ' -f2 ${fsdg_file} | sort
            echo ${output_path}/check2_temp.txt
            echo ${output_path}/check1_temp.txt

            if diff ${output_path}/check1_temp.txt ${output_path}/check2_temp.txt > /dev/null; then
                echo "Files are identical. Continuing..."
                #rm ${output_path}/check2_temp.txt
                #rm ${output_path}/check1_temp.txt
            else
                echo "################################################"
                echo "Files differ! Aborting script."
                echo "################################################"
                #rm ${output_path}/check2_temp.txt
                #rm ${output_path}/check1_temp.txt
                exit 1
            fi
                
            if [[ "$fsdg_file" == *"average"* ]]; then  # forthe group level do one type of model 
                if [[ "$do_fitting" == "1" ]]; then
                        mri_glmfit \
                            --y ${output_path}/${hemi}.long.thickness-${long_parameter}.fwhm${fwhm_value}.fsaverage.4d.mgh  \
                            --fsgd ${fsdg_file} \
                            --C ${output_path}/average.mtx \
                            --C ${output_path}/age.mtx \
                            --C ${output_path}/sex.mtx \
                            --surf fsaverage ${hemi} \
                            --cortex --eres-save \
                            --glmdir ${output_path}/glm_average_group_${hemi} --debug #--checkopts 

                    if [ $correction_with_MCZ == "1" ]; then
                        mri_glmfit-sim --glmdir ${output_path}/glm_average_group_${hemi} \
                                    --y  ${output_path}/${hemi}.long.thickness-${long_parameter}.fwhm${fwhm_value}.fsaverage.4d.mgh \
                                    --cwp 0.05 \
                                    --2spaces \
                                    --mczsim ${cluster_forming} abs \
                                    --overwrite 
                    elif [ $correction_with_permutation == "1" ]; then
                        # correction for multiple comparison, usign both thinning (pos) and tickening (neg)
                        mri_glmfit-sim --glmdir ${output_path}/glm_average_group_${hemi} \
                                        --y  ${output_path}/${hemi}.long.thickness-${long_parameter}.fwhm${fwhm_value}.fsaverage.4d.mgh \
                                        --cwp 0.05 \
                                        --perm ${NPERM} ${TH} abs \
                                        --perm-resid \
                                        --bg ${nthreads} \
                                        --2spaces \
                                        --overwrite #--seed 42 if you want reproducibility you can set the seed but then you must use a single thread otherwise it gives problem but of course is way slower

                       # mri_glmfit-sim --glmdir ${output_path}/glm_average_group_${hemi} \
                       #                 --y  ${output_path}/${hemi}.long.thickness-${long_parameter}.fwhm${fwhm_value}.fsaverage.4d.mgh \
                       #                 --perm ${NPERM} ${TH} neg \
                       #                 --perm-resid \
                       #                 --2spaces \
                       #                 --bg ${nthreads}  \
                       #                 --overwrite 
                    fi 
                   
                fi 
            elif [[ "$fsdg_file" == *"fenotipo"* ]] || [[ "$fsdg_file" == *"MCI"* ]]; then 
                if [[ "$do_fitting" == "1" ]]; then
                    mri_glmfit \
                        --y ${output_path}/${hemi}.long.thickness-${long_parameter}.fwhm${fwhm_value}.fsaverage.4d.mgh  \
                        --fsgd ${fsdg_file} \
                        --C ${output_path}/between_group_difference.mtx \
                        --C ${output_path}/interaction_group_age.mtx \
                        --C ${output_path}/interaction_group_sex.mtx \
                        --C ${output_path}/interaction_group_age_sex.mtx \
                        --surf fsaverage ${hemi} \
                        --cortex --eres-save \
                        --glmdir ${output_path}/glm_average_group_${hemi} --debug
                    if [ $correction_with_MCZ == "1" ]; then

                        mri_glmfit-sim --glmdir ${output_path}/glm_average_group_${hemi} \
                                    --y ${output_path}/${hemi}.long.thickness-${long_parameter}.fwhm${fwhm_value}.fsaverage.4d.mgh \
                                    --cwp 0.05 \
                                    --2spaces \
                                    --mczsim ${cluster_forming} abs \
                                    --overwrite 

                    elif [ $correction_with_permutation == "1" ]; then
                    
                        mri_glmfit-sim --glmdir ${output_path}/glm_average_group_${hemi} \
                                        --y  ${output_path}/${hemi}.long.thickness-${long_parameter}.fwhm${fwhm_value}.fsaverage.4d.mgh \
                                        --perm ${NPERM} ${TH} abs \
                                        --perm-resid \
                                        --bg ${nthreads} \
                                        --2spaces \
                                        --overwrite #--seed 42 if you want reproducibility you can set the seed but then you must use a single thread otherwise it gives problem but of course is way slower

                        #mri_glmfit-sim --glmdir ${output_path}/glm_average_group_${hemi} \
                        #                --y  ${output_path}/${hemi}.long.thickness-${long_parameter}.fwhm${fwhm_value}.fsaverage.4d.mgh \
                        #                --perm ${NPERM} ${TH} neg \
                        #                --perm-resid \
                        #                --2spaces \
                        #                --bg ${nthreads} \
                        #                --overwrite 

                    fi 
                fi 
            fi 
        done    
        # group average, first plot is uncorrected with threhsold low_p (i.e. 0.01) and the second plot is the corrected with permutation testing
        # I am showing in the permutatin testing only the negative becasue the hypothesis is on the cortical thinning (i.e. assottigliamento)

        if [[ "$do_plotting" == "1" ]]; then 

            if  [ $correction_with_permutation == "1" ]; then
                if [[ "$fsdg_file" == *"average"* ]]; then 
                    for analysis in average; do #{average,age,protocol}; do 
                        echo "uncorrected $analysis for $fsdg_file"
                        freeview -f $SUBJECTS_DIR/fsaverage/surf/lh.inflated:overlay=${output_path}/glm_average_group_lh/${analysis}/sig.mgh:overlay_threshold=${low_p},${high_p}:annot=$SUBJECTS_DIR/fsaverage/label/lh.aparc.annot:annot_outline=1 \
                                    $SUBJECTS_DIR/fsaverage/surf/rh.inflated:overlay=${output_path}/glm_average_group_rh/${analysis}/sig.mgh:overlay_threshold=${low_p},${high_p}:annot=$SUBJECTS_DIR/fsaverage/label/rh.aparc.annot:annot_outline=1
                    done 

                    for analysis in average; do #,age,protocol}; do 
                        for sign in abs; do #{pos,neg}; do
                            echo "corrected $analysis $sign for $output_stats-based correction for $fsdg_file" 
                            #less ${output_path}/glm_average_group_lh/${analysis}/perm.th13.${sign}.sig.${output_stats}.summary
                            mris_calc ${output_path}/glm_average_group_lh/${analysis}/perm.th13.${sign}.sig.${output_stats}.mgh max 
                            mris_calc ${output_path}/glm_average_group_lh/${analysis}/perm.th13.${sign}.sig.${output_stats}.mgh min
                            mris_calc ${output_path}/glm_average_group_rh/${analysis}/perm.th13.${sign}.sig.${output_stats}.mgh max
                            mris_calc ${output_path}/glm_average_group_rh/${analysis}/perm.th13.${sign}.sig.${output_stats}.mgh min

                            freeview -f $SUBJECTS_DIR/fsaverage/surf/lh.pial:overlay=${output_path}/glm_average_group_lh/${analysis}/perm.th13.${sign}.sig.${output_stats}.mgh:overlay_threshold=${low_p_corrected},${high_p_corrected}:annot=$SUBJECTS_DIR/fsaverage/label/lh.aparc.annot:annot_outline=1 \
                                $SUBJECTS_DIR/fsaverage/surf/rh.pial:overlay=${output_path}/glm_average_group_rh/${analysis}/perm.th13.${sign}.sig.${output_stats}.mgh:overlay_threshold=${low_p_corrected},${high_p_corrected}:annot=$SUBJECTS_DIR/fsaverage/label/rh.aparc.annot:annot_outline=1
                        done 
                    done 

                elif [[ "$fsdg_file" == *"fenotipo"* ]] || [[ "$fsdg_file" == *"MCI"* ]]; then  

                    for analysis in between_group_difference; do # ,interaction_group_age,interaction_group_protocol,interaction_group_age_protocol}; do 
                        echo "uncorrected $analysis for $fsdg_file"
                        freeview -f $SUBJECTS_DIR/fsaverage/surf/lh.pial:overlay=${output_path}/glm_average_group_lh/${analysis}/sig.mgh:overlay_threshold=${low_p},${high_p}:annot=$SUBJECTS_DIR/fsaverage/label/lh.aparc.annot:annot_outline=1 \
                                    $SUBJECTS_DIR/fsaverage/surf/rh.pial:overlay=${output_path}/glm_average_group_rh/${analysis}/sig.mgh:overlay_threshold=${low_p},${high_p}:annot=$SUBJECTS_DIR/fsaverage/label/rh.aparc.annot:annot_outline=1
                    done 

                    for analysis in between_group_difference; do #,interaction_group_age,interaction_group_protocol,interaction_group_age_protocol}; do 
                        for sign in abs; do #{pos,neg}; do
                            echo "corrected $analysis $sign for $output_stats-based correction for $fsdg_file" 
                            #less ${output_path}/glm_average_group_lh/${analysis}/perm.th13.${sign}.sig.${output_stats}.superm.th13.${sign}.sig.${output_stats}.mghmmary
                            mris_calc ${output_path}/glm_average_group_lh/${analysis}/perm.th13.${sign}.sig.${output_stats}.mgh max 
                            mris_calc ${output_path}/glm_average_group_lh/${analysis}/perm.th13.${sign}.sig.${output_stats}.mgh min
                            mris_calc ${output_path}/glm_average_group_rh/${analysis}/perm.th13.${sign}.sig.${output_stats}.mgh max
                            mris_calc ${output_path}/glm_average_group_rh/${analysis}/perm.th13.${sign}.sig.${output_stats}.mgh min

                            freeview -f $SUBJECTS_DIR/fsaverage/surf/lh.pial:overlay=${output_path}/glm_average_group_lh/${analysis}/perm.th13.${sign}.sig.${output_stats}.mgh:overlay_threshold=${low_p_corrected},${high_p_corrected}:annot=$SUBJECTS_DIR/fsaverage/label/lh.aparc.annot:annot_outline=1 \
                                        $SUBJECTS_DIR/fsaverage/surf/rh.pial:overlay=${output_path}/glm_average_group_rh/${analysis}/perm.th13.${sign}.sig.${output_stats}.mgh:overlay_threshold=${low_p_corrected},${high_p_corrected}:annot=$SUBJECTS_DIR/fsaverage/label/rh.aparc.annot:annot_outline=1 \
                                        --viewport 3d
                        done 
                    done

                fi 
            elif [ $correction_with_MCZ == "1" ]; then 
                if [[ "$fsdg_file" == *"average"* ]]; then 
                    for analysis in average; do #{average,age,protocol}; do 
                        echo "uncorrected $analysis for $fsdg_file"
                        freeview -f $SUBJECTS_DIR/fsaverage/surf/lh.pial:overlay=${output_path}/glm_average_group_lh/${analysis}/sig.mgh:overlay_threshold=${low_p},${high_p}:annot=$SUBJECTS_DIR/fsaverage/label/lh.aparc.annot:annot_outline=1 \
                                    $SUBJECTS_DIR/fsaverage/surf/rh.pial:ove
                    done 

                    for analysis in average; do #,age,protocol}; do 
                        for sign in abs; do
                            echo "corrected $analysis $sign for $output_stats-based correction for $fsdg_file" 
                            #less ${output_path}/glm_average_group_lh/${analysis}/perm.th13.${sign}.sig.${output_stats}.summary
                            echo ${output_path}/glm_average_group_lh/${analysis}/cache.th30.abs.sig.cluster.summary
                            freeview -f $SUBJECTS_DIR/fsaverage/surf/lh.pial:overlay=${output_path}/glm_average_group_lh/${analysis}/cache.th30.abs.sig.${output_stats}.mgh:overlay_threshold=${low_p_corrected},${high_p_corrected}:annot=$SUBJECTS_DIR/fsaverage/label/lh.aparc.annot:annot_outline=1 \
                                $SUBJECTS_DIR/fsaverage/surf/rh.pial:overlay=${output_path}/glm_average_group_rh/${analysis}/cache.th30.abs.sig.${output_stats}.mgh:overlay_threshold=${low_p_corrected},${high_p_corrected}:annot=$SUBJECTS_DIR/fsaverage/label/rh.aparc.annot:annot_outline=1
                        done 
                    done 

                elif [[ "$fsdg_file" == *"fenotipo"* ]] || [[ "$fsdg_file" == *"MCI"* ]]; then  

                    for analysis in between_group_difference; do # ,interaction_group_age,interaction_group_protocol,interaction_group_age_protocol}; do 
                        echo "uncorrected $analysis for $fsdg_file"
                        freeview  -f $SUBJECTS_DIR/fsaverage/surf/lh.pial:overlay=${output_path}/glm_average_group_lh/${analysis}/sig.mgh:overlay_threshold=${low_p},${high_p}:annot=$SUBJECTS_DIR/fsaverage/label/lh.aparc.annot:annot_outline=1 \
                                     $SUBJECTS_DIR/fsaverage/surf/rh.pial:overlay=${output_path}/glm_average_group_rh/${analysis}/sig.mgh:overlay_threshold=${low_p},${high_p}:annot=$SUBJECTS_DIR/fsaverage/label/rh.aparc.annot:annot_outline=1
                    done 

                    for analysis in between_group_difference; do #,interaction_group_age,interaction_group_protocol,interaction_group_age_protocol}; do 
                        for sign in abs; do
                            echo "corrected $analysis $sign for $output_stats-based correction for $fsdg_file" 
                            cat ${output_path}/glm_average_group_lh/${analysis}/cache.th30.abs.sig.cluster.summary
                            freeview -f $SUBJECTS_DIR/fsaverage/surf/lh.pial:overlay=${output_path}/glm_average_group_lh/${analysis}/cache.th30.abs.sig.${output_stats}.mgh:overlay_threshold=${low_p_corrected},${high_p_corrected}:annot=$SUBJECTS_DIR/fsaverage/label/lh.aparc.annot:annot_outline=1 \
                                        $SUBJECTS_DIR/fsaverage/surf/rh.pial:overlay=${output_path}/glm_average_group_rh/${analysis}/cache.th30.abs.sig.${output_stats}.mgh:overlay_threshold=${low_p_corrected},${high_p_corrected}:annot=$SUBJECTS_DIR/fsaverage/label/rh.aparc.annot:annot_outline=1 \
                                        --viewport 3d
                        done 
                    done

                fi 
            fi 
        fi 
    done
fi
    
    #if [ $fdr_correction == "0" ]; then 
        # take vertex that are signficant at p < 0.05 i.e. -log10(p) > 1.3 [rememebr 1.3 0.05, 2.0 0.01, 3.0 0.001, 3.7 0.0002]
    #    mri_binarize --i ${output_path}/glm_group_analysis_lh/lh_group/sig.mgh --min 3.0 --binval 1 --o ${output_path}/glm_group_analysis_lh/lh_group/sig_masked.mgh
    #    mri_binarize --i ${output_path}/glm_group_analysis_rh/rh_group/sig.mgh --min 3.0 --binval 1 --o ${output_path}/glm_group_analysis_rh/rh_group/sig_masked.mgh

    #    mri_mask ${output_path}/glm_group_analysis_lh/lh_group/gamma.mgh ${output_path}/glm_group_analysis_lh/lh_group/sig_masked.mgh ${output_path}/glm_group_analysis_lh/lh_group/gamma_masked.mgh
    #    mri_mask ${output_path}/glm_group_analysis_rh/rh_group/gamma.mgh ${output_path}/glm_group_analysis_rh/rh_group/sig_masked.mgh  ${output_path}/glm_group_analysis_rh/rh_group/gamma_masked.mgh
                
        #mris_calc ${output_path}/glm_group_analysis_rh/rh_group/gamma_masked.mgh max
        #mris_calc ${output_path}/glm_group_analysis_lh/lh_group/gamma_masked.mgh max
        #mris_calc ${output_path}/glm_group_analysis_rh/rh_group/gamma_masked.mgh min 
        #mris_calc ${output_path}/glm_group_analysis_lh/lh_group/gamma_masked.mgh min

        # p = 0.01 is 2, p < 0.05 is 1.3
        #freeview -f $SUBJECTS_DIR/fsaverage/surf/lh.pial:overlay=${output_path}/glm_group_analysis_lh/lh_group/sig_masked.mgh:overlay_threshold=0.5,1 \
        #         -f $SUBJECTS_DIR/fsaverage/surf/lh.pial:overlay=${output_path}/glm_group_analysis_lh/lh_group/sig.mgh \
        #         -f $SUBJECTS_DIR/fsaverage/surf/rh.pial:overlay=${output_path}/glm_group_analysis_rh/rh_group/sig_masked.mgh:overlay_threshold=0.5,1 \
        #         -f $SUBJECTS_DIR/fsaverage/surf/rh.pial:overlay=${output_path}/glm_group_analysis_rh/rh_group/sig.mgh  --viewport 3d #--colorscale 
        #freeview -f $SUBJECTS_DIR/fsaverage/surf/lh.pial:overlay=${output_path}/glm_group_analysis_lh/lh_group/gamma_masked.mgh:overlay_threshold=0.001,3.5 \
        #         -f $SUBJECTS_DIR/fsaverage/surf/rh.pial:overlay=${output_path}/glm_group_analysis_rh/rh_group/gamma_masked.mgh:overlay_threshold=0.001,3.5 \
        #          --viewport 3d --colorscale 

  