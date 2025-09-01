import pandas as pd
import argparse 
import sys
from loguru import logger
from utils.utils import load_img, save_img, parse_yaml, str_to_bool, LoguruStreamHandler
import numpy as np
import os
from unidecode import unidecode
from scipy.stats import shapiro
from scipy.stats import mannwhitneyu
from scipy.stats import ttest_ind
from scipy.stats import chi2_contingency
from scipy.stats import fisher_exact
from termcolor import colored, cprint

def parse_args(args):
    """parse argument"""
    parser = argparse.ArgumentParser(epilog = "This script is done to define the file for performing statistics on the Salerno MSA study. Longitudinal freesurer run")
    parser.add_argument("-pf", "--patient_file", required=True, type = str, help = "This file is outputted from the dicom2bids.py script. This is the file with the patient information")
    parser.add_argument("-doqdec", "--do_qdecfile", required=True, type = str_to_bool, default = False, help = "if true prepares and output in the stats results foldet the qdec file to run longitudinal mris_long_slope \
                                                                                                                 rate of change thickness")
    parser.add_argument("-dofsdg", "--do_fsdgfile", required=True, type = str_to_bool, default = False, help = "if True prepares and output the fsdg file for different design: groupe level analysis and between group analysis")
    parser.add_argument("-protocol", "--protocol", required=True, type = str, help = "different way to treat the protocol variable when defining the GLM." \
                                                                                    "'all' use all three protocols, 'two' remove the subjects from the " \
                                                                                    "smallest protocol group N = 4 and 'merge' the smallest protocol with the protcol with N = 21 MSA that is the")

    return parser.parse_args()
 
def main(argv):
    try:
        print('this is the main')
        args = parse_args(argv)
        PAZ_data = pd.read_excel(args.patient_file,sheet_name="Sheet1",index_col=None, engine='openpyxl')    
        PAZ_data.columns = [unidecode(col) for col in PAZ_data.columns]
        PAZ_data['UMSARSIdelta'] = PAZ_data['UMSARS I.1'] - PAZ_data['UMSARS I']
        PAZ_data['UMSARSIIdelta'] = PAZ_data['UMSARS II.1'] - PAZ_data['UMSARS II']
        print(PAZ_data.columns)

        if args.do_qdecfile == True:
            for number_oftimepoints in [2,3,4]:  # this is the loop for the patients that have longitudinal runs
                #print(number_oftimepoints)
                # Step 1: Group and collect unique sessions per subject
                if number_oftimepoints == 2:

                        session_counts = PAZ_data.groupby('Subj_ID')['Session'].apply(set)
                        subjects_with_T0_T1 = session_counts[session_counts.apply(lambda x: 'ses-T0' in x and 'ses-T1' in x)].index 
                        print("################")
                        #print(subjects_with_T0_T1)
                        filtered_df = PAZ_data[PAZ_data['Subj_ID'].isin(subjects_with_T0_T1)]
                        filtered_df = filtered_df[~filtered_df['Session'].str.contains('ses-T2')]
                        filtered_df = filtered_df[~filtered_df['Session'].str.contains('ses-T3')]
                        #print(filtered_df)
                        output_filename="two"

                elif  number_oftimepoints == 3:

                        session_counts = PAZ_data.groupby('Subj_ID')['Session'].apply(set)
                        subjects_with_T0_T1_T2 = session_counts[session_counts.apply(lambda x: 'ses-T0' in x and 'ses-T1' and 'ses-T2' in x)].index 
                        filtered_df = PAZ_data[PAZ_data['Subj_ID'].isin(subjects_with_T0_T1_T2)]
                        filtered_df = filtered_df[~filtered_df['Session'].str.contains('ses-T3')]
                        output_filename="three" 

                elif  number_oftimepoints == 4:
                    
                        session_counts = PAZ_data.groupby('Subj_ID')['Session'].apply(set)
                        subjects_with_T0_T1_T2_T3 = session_counts[session_counts.apply(lambda x: 'ses-T0' in x and 'ses-T1' and 'ses-T2' and 'ses-T3' in x)].index 
                        filtered_df = PAZ_data[PAZ_data['Subj_ID'].isin(subjects_with_T0_T1_T2_T3)]
                        output_filename="four"          

                # Step 2 populate the new dataframe 
                qdec = pd.DataFrame(columns=['fsid', 'fsid-base', 'years', 'Eta', 'sex','diff']) # then you drop diff column 
                qdec['sex'] = filtered_df['Sex'].map({'M': 0, 'F': 1})
                qdec['diff'] = filtered_df['T1_T0_Diff_Days']/365
                qdec['fsid-base'] = filtered_df['Subj_ID'] + '_'
                qdec['fsid'] = filtered_df['Subj_ID'] + '_' + filtered_df['Session']
                qdec['years'] = round(filtered_df['years']/365,3)
                qdec['age_old'] = filtered_df['Eta'] #.str.rstrip('Y'), errors='coerce').astype('Int64')

                #adjust the age value
                qdec['subject'] = qdec['fsid'].str.extract(r'(sub-\d+)_')
                t0_values = qdec[qdec['fsid'].str.contains("ses-T0")].copy()
                t0_values['age_plus_diff'] = t0_values['age_old'] + t0_values['diff']
                t0_map = t0_values.set_index('subject')['age_plus_diff']
                qdec['age_T0_plus_diff'] = qdec['subject'].map(t0_map)
                #print(qdec)
                qdec['age'] = qdec.apply(lambda row: row['age_T0_plus_diff'] if "ses-T1" in row['fsid'] else row['age_old'],axis=1)
                        
                #print(qdec)
                #print(qdec.shape)
                qdec.drop(columns=['age','sex','diff','age_old','subject','age_T0_plus_diff'], axis = 1, inplace = True)
                qdec.to_csv(f"/home/riccardo/codici_progetti/Salerno/qdec/{output_filename}_timepoints_allMSA_long.qdec.table.dat", sep=" ", index=False)

            ######## this is the part for the helthy controls which instead are cross-sectional. do this when they will be finished cross sectionally (learn how to use this in the cross sectional setting paz vs controls)
        if args.do_fsdgfile == True:
            print(PAZ_data.shape)

            # getting only subjects with two timepoints
            session_counts = PAZ_data.groupby('Subj_ID')['Session'].apply(set)
            subjects_with_T0_T1 = session_counts[session_counts.apply(lambda x: 'ses-T0' in x and 'ses-T1' in x)].index 
            print(len(subjects_with_T0_T1))
            PAZ_data = PAZ_data[PAZ_data['Subj_ID'].isin(subjects_with_T0_T1)]
            PAZ_data = PAZ_data[~PAZ_data['Session'].str.contains('ses-T2')]
            PAZ_data = PAZ_data[~PAZ_data['Session'].str.contains('ses-T3')]

            # remove subject with 1.5 T and the outliers
            PAZ_data = PAZ_data[~PAZ_data["Subj_ID"].str.contains("sub-18", na=False)] # 1.5 T at baseline (Philips)
            PAZ_data = PAZ_data[~PAZ_data["Subj_ID"].str.contains("sub-45", na=False)] # outlier
            PAZ_data = PAZ_data[~PAZ_data["Subj_ID"].str.contains("sub-53", na=False)] # outlier

            # add a categorical variables for the protocol
            unique_times = sorted(PAZ_data['RepetitionTime'].dropna().unique())

            # Assign a group name to each unique RepetitionTime
            group_map = {time: f'Group {chr(65 + i)}' for i, time in enumerate(unique_times)}
            PAZ_data['protocol'] = PAZ_data['RepetitionTime'].map(group_map)

            # maintain only the data important for the fsdg file
            PAZ_data  = PAZ_data[['Subj_ID','Session','Fenotipo','MCI_T0','Eta','Sex','protocol','UMSARSIdelta','UMSARSIIdelta', \
                                   'T1_T0_Diff_Days','Durata malattia (anni)','data esame BL','data esame FUP 1','years', \
                                    'UMSARS I','UMSARS I.1','UMSARS II','UMSARS II.1','UMSARS IV','UMSARS IV.1']] # reordering 
            PAZ_data['Sex'] = PAZ_data['Sex'].map({'M': 0, 'F': 1})
            PAZ_data = PAZ_data[~PAZ_data['Session'].str.contains('ses-T1')]

            # convert the protocol A B and C to 1 2 3
            mapping = {"Group A": 1, "Group B": 2, "Group C": 3}
            PAZ_data["protocol"] = PAZ_data["protocol"].map(mapping)
            if args.protocol == 'all':
                out_name = 'allprotocols_nosubj18nosubj45nosubj53'  # subj 18 is the one with 1.5 T the others two are outliers whih may skwew the analysis
            elif args.protocol == 'two':
                out_name = 'twoprotocols_nosubj18nosubj45nosubj53'
                PAZ_data = PAZ_data[PAZ_data['protocol'] != 1] # 1 is the label for the protocol group with less people
            elif args.protocol == 'merge':
                out_name = 'mergeprotocols_nosubj18nosubj45nosubj53'
                PAZ_data['protocol'].replace(1,2, inplace=True)
                
            PAZ_data["MCI_T0"] = PAZ_data["MCI_T0"].astype(int)
            PAZ_data['Age'] = PAZ_data['Eta'] #.str.rstrip('Y'), errors='coerce').astype('Int64')

            # final edits to the dataframe
            PAZ_data.drop(columns=['Session'], axis = 1, inplace = True)
            PAZ_data['Subj_ID'] = PAZ_data['Subj_ID'] + '_'
            #print(PAZ_data)
            ###############################################
            # prepare fsdg for group level with covariates 
            ###############################################
            out_dir=f"./stats_results/MSA_group_average_{out_name}"
            out_fsdg_name="MSA_average_group_age_sex.fsdg"
            if os.path.isdir(out_dir) == False:
                os.makedirs(out_dir)

            with open(f"{out_dir}/{out_fsdg_name}", "w") as f:
                f.write("GroupDescriptorFile 1\n")
                f.write("Title MSA CHANGE\n")
                f.write("Class MSA\n")
                f.write("Variables Age Sex\n")

                for _, row in PAZ_data.iterrows():
                    line = f"Input {row['Subj_ID']} MSA {row['Age']} {row['Sex']}\n"
                    f.write(line)

            # contrast files for the group average
            np.savetxt(f"{out_dir}/average.mtx", np.array([[1, 0, 0,]]), fmt="%.1f")
            np.savetxt(f"{out_dir}/age.mtx", np.array([[0, 1, 0]]), fmt="%.1f")
            np.savetxt(f"{out_dir}/sex.mtx", np.array([[0, 0, 1]]), fmt="%.1f")
            #np.savetxt(f"{out_dir}/protocol.mtx", np.array([[0, 0, 1]]), fmt="%.1f")

            ###############################################
            # prepare fsdg for between group comparison MSAp vs MSAc
            ###############################################
            out_dir=f"./stats_results/MSA_between_group_fenotipo_{out_name}"
            out_fsdg_name="MSA_between_group_fenotipo.fsdg"
            if os.path.isdir(out_dir) == False:
                os.makedirs(out_dir)

            with open(f"{out_dir}/{out_fsdg_name}", "w") as f:
                f.write("GroupDescriptorFile 1\n")
                f.write("Title MSA CHANGE MSAp vs MSAc\n")
                f.write("Class MSAP\n")
                f.write("Class MSAC\n")
                f.write("Variables Age Sex\n")

                for _, row in PAZ_data.iterrows():
                    if row['Fenotipo'] == 'C':
                        line = f"Input {row['Subj_ID']} MSAC {row['Age']} {row['Sex']}\n"
                    elif row['Fenotipo'] == 'P':
                        line = f"Input {row['Subj_ID']} MSAP {row['Age']} {row['Sex']}\n"
                    f.write(line)

            # contrast files for the group average
            np.savetxt(f"{out_dir}/between_group_difference.mtx", np.array([[1, -1, 0, 0, 0, 0]]), fmt="%.1f") # this is group 1 > group2 i.e. MSAp > MSAc
            np.savetxt(f"{out_dir}/interaction_group_age.mtx", np.array([[0, 0, 1, -1, 0, 0]]), fmt="%.1f") # interaction between group and age
            np.savetxt(f"{out_dir}/interaction_group_sex.mtx", np.array([[0, 0, 0, 0, 1, -1]]), fmt="%.1f") # interaction between group and protocol
            np.savetxt(f"{out_dir}/interaction_group_age_sex.mtx", np.array([[0, 0, 1, -1, 0, 0], [0, 0, 0, 0, 1, -1]]), fmt="%.1f") # interaction between group and age and protocol

            ###############################################
            # prepare fsdg for between group comparison MCIyes vs MCIno
            ###############################################
            out_dir=f"./stats_results/MSA_between_group_MCI_{out_name}"
            out_fsdg_name="MSA_between_group_MCI.fsdg"
            if os.path.isdir(out_dir) == False:
                os.makedirs(out_dir)
            
            with open(f"{out_dir}/{out_fsdg_name}", "w") as f:
                f.write("GroupDescriptorFile 1\n")
                f.write("Title MSA CHANGE MCIyes vs MCIno\n")
                f.write("Class MCIyes\n")
                f.write("Class MCIno\n")
                f.write("Variables Age Sex\n")

                for _, row in PAZ_data.iterrows():
                    if row['MCI_T0'] == 0:
                        line = f"Input {row['Subj_ID']} MCIno {row['Age']} {row['Sex']}\n"
                    elif row['MCI_T0'] == 1:
                        line = f"Input {row['Subj_ID']} MCIyes {row['Age']} {row['Sex']}\n"

                    f.write(line)

            # contrast files for the group average
            np.savetxt(f"{out_dir}/between_group_difference.mtx", np.array([[1, -1, 0, 0, 0, 0]]), fmt="%.1f") # this is group 1 > group2 i.e. MSAp > MSAc
            np.savetxt(f"{out_dir}/interaction_group_age.mtx", np.array([[0, 0, 1, -1, 0, 0]]), fmt="%.1f") # interaction between group and age
            np.savetxt(f"{out_dir}/interaction_group_sex.mtx", np.array([[0, 0, 0, 0, 1, -1]]), fmt="%.1f") # interaction between group and protocol
            np.savetxt(f"{out_dir}/interaction_group_age_sex.mtx", np.array([[0, 0, 1, -1, 0, 0], [0, 0, 0, 0, 1, -1]]), fmt="%.1f") # interaction between group and age and protocol

            ###############################################
            # prepare fsdg for between group comparison MCIyes vs MCIno in the group of MSAc only
            ###############################################
            PAZ_data_MSAc = PAZ_data[PAZ_data['Fenotipo'].isin(['C'])]
            out_dir=f"./stats_results/MSAcerebellar_between_group_MCI_{out_name}"
            out_fsdg_name="MSAcerebellar_between_group_MCI.fsdg"
            if os.path.isdir(out_dir) == False:
                os.makedirs(out_dir)
            with open(f"{out_dir}/{out_fsdg_name}", "w") as f:
                f.write("GroupDescriptorFile 1\n")
                f.write("Title MSAcerebellar CHANGE MCIyes vs MCIno\n")
                f.write("Class MCIyes\n")
                f.write("Class MCIno\n")
                f.write("Variables Age Sex\n")

                for _, row in PAZ_data_MSAc.iterrows():  # user here the filtered dataset with only the MSA cerebellar phenotype :)
                    if row['MCI_T0'] == 0:
                        line = f"Input {row['Subj_ID']} MCIno {row['Age']} {row['Sex']}\n"
                    elif row['MCI_T0'] == 1:
                        line = f"Input {row['Subj_ID']} MCIyes {row['Age']} {row['Sex']}\n"

                    f.write(line)

            #print(PAZ_data)
            #print("########")
            #print(PAZ_data_MSAc)
            # contrast files for the group average
            np.savetxt(f"{out_dir}/between_group_difference.mtx", np.array([[1, -1, 0, 0, 0, 0]]), fmt="%.1f") # this is group 1 > group2 i.e. MSAp > MSAc
            np.savetxt(f"{out_dir}/interaction_group_age.mtx", np.array([[0, 0, 1, -1, 0, 0]]), fmt="%.1f") # interaction between group and age
            np.savetxt(f"{out_dir}/interaction_group_sex.mtx", np.array([[0, 0, 0, 0, 1, -1]]), fmt="%.1f") # interaction between group and protocol
            np.savetxt(f"{out_dir}/interaction_group_age_sex.mtx", np.array([[0, 0, 1, -1, 0, 0], [0, 0, 0, 0, 1, -1]]), fmt="%.1f") # interaction between group and age and protocol

            ###############################################
            # prepare fsdg for within group comparison MSAc and MSAp with the spc of subcortical volumetric change 
            ###############################################

            for hemi in ['lh','rh']:
                for feno in ['C','P','all']:
                    if 'C' in feno or 'P' in feno:
                        PAZ_feno = PAZ_data[PAZ_data['Fenotipo'].isin([feno])] # otherwise on all subjects 
                    else:
                        PAZ_feno = PAZ_data

                    for subcortical in ['Putamen','brain_stem','Pallidum','CerebellumCortex','CerebellumWhiteMatter','CerebellumTot','Whole_brainstem','SCP','Pons','BrainstemPons']:
                        if "Whole_brainstem" in subcortical or "SCP" in subcortical or "Pons" in subcortical:  # load brainstem data whihc are done with the custom sgmentation of brainstem rather than the aseg of freesurfer
                            if "BrainstemPons" in subcortical:
                                vols = pd.read_csv(f"/home/riccardo/codici_progetti/Salerno/stats_results/correlations/spc_{hemi}_brainstem_MSA{feno}_{subcortical}.csv", header=None, names=["Pons", "Whole_brainstem", "Subject"])
                            else:
                                vols = pd.read_csv(f"/home/riccardo/codici_progetti/Salerno/stats_results/correlations/spc_{hemi}_brainstem_MSA{feno}_{subcortical}.csv", header=None, names=[f"{subcortical}", "Subject"])
                        else:
                            vols = pd.read_csv(f"/home/riccardo/codici_progetti/Salerno/stats_results/correlations/spc_{hemi}_MSA{feno}_{subcortical}.csv", header=None, names=[f"{subcortical}", "Subject"])

                        vols['Subject'] = vols['Subject'].str[:7]
                        vols.rename(columns={'Subject':'Subj_ID'}, inplace = True)

                        col1 = vols['Subj_ID'].reset_index(drop=True)
                        col2 = PAZ_feno['Subj_ID'].reset_index(drop=True)
                        are_equal = col1.equals(col2)
                        if are_equal == False:
                            print('files are not matching')
                            sys.exit()
                        #else: print('everything ok')
                        
                        merged = pd.merge(PAZ_feno, vols, on = 'Subj_ID', how = 'inner')

                        if "BrainstemPons" in subcortical and "all" in feno:
                            print("################")
                            print(merged)
                            
                        out_dir=f"./stats_results/correlations/MSA{feno}_subcort_correlation_{subcortical}_{out_name}"
                        out_fsdg_name=f"MSA{feno}_subcort_correlation_{hemi}_{subcortical}.fsdg"

                        if os.path.isdir(out_dir) == False:
                            os.makedirs(out_dir)
                        with open(f"{out_dir}/{out_fsdg_name}", "w") as f:
                            f.write("GroupDescriptorFile 1\n")
                            f.write(f"Title MSA{feno} correlation subcortical\n")
                            f.write(f"Class MSA{feno}\n")
                            if "BrainstemPons" in subcortical:
                                f.write(f"Variables Age Sex {hemi}_Pons {hemi}_Whole_brainstem\n")
                            else:
                                f.write(f"Variables Age Sex {hemi}_{subcortical}\n")

                            for _, row in merged.iterrows():  # user here the filtered dataset with only the MSA cerebellar phenotype :)
                                #print(row[f'{subcortical}'])
                                if "BrainstemPons" in subcortical:
                                    print(row)
                                    line = f"Input {row['Subj_ID']} MSA{feno} {row['Age']} {row['Sex']} {row['Pons']} {row['Whole_brainstem']}\n"
                                else:
                                    line = f"Input {row['Subj_ID']} MSA{feno} {row['Age']} {row['Sex']} {row[f'{subcortical}']}\n"
                                f.write(line)

                        if 'BrainstemPons' in subcortical:
                            np.savetxt(f"{out_dir}/correlation_subcortical_{hemi}_{subcortical}_Pons.mtx", np.array([[0, 0, 0, 1, 0]]), fmt="%.1f") 
                            np.savetxt(f"{out_dir}/correlation_subcortical_{hemi}_{subcortical}_Wholebrainstem.mtx", np.array([[0, 0, 0, 0, 1]]), fmt="%.1f") 
                        else:
                            np.savetxt(f"{out_dir}/correlation_subcortical_{hemi}_{subcortical}.mtx", np.array([[0, 0, 0, 1]]), fmt="%.1f") 

            ###############################################
            # prepare fsdg for within group comparison MSAc, MSAp, MSAall correlated with UMSARSvariable
            ###############################################
            for hemi in ['lh','rh']:
                for feno in ['C','P','all']:
                    if 'C' in feno or 'P' in feno:
                        PAZ_feno = PAZ_data[PAZ_data['Fenotipo'].isin([feno])] # otherwise on all subjects 
                    else:
                        PAZ_feno = PAZ_data

                    for scale in ['UMSARSIdelta','UMSARSIIdelta']: 
                        out_dir=f"./stats_results/correlations/MSA{feno}_correlation_with_{scale}_{out_name}" # 
                        out_fsdg_name=f"MSA{feno}_correlation_scale_{scale}.fsdg"
                        if os.path.isdir(out_dir) == False:
                            print(feno)
                            os.makedirs(out_dir)
                        with open(f"{out_dir}/{out_fsdg_name}", "w") as f:
                            f.write("GroupDescriptorFile 1\n")
                            f.write(f"Title MSA{feno} correaltion {scale}\n")
                            f.write(f"Class MSA{feno}\n")
                            f.write(f"Variables Age Sex {scale}\n")

                            for _, row in PAZ_feno.iterrows():  # user here the filtered dataset with only the MSA cerebellar phenotype :)
                                line = f"Input {row['Subj_ID']} MSA{feno} {row['Age']} {row['Sex']} {row[scale]}\n"
                                f.write(line)
                    
                            np.savetxt(f"{out_dir}/correlation_subcortical_{hemi}_{scale}.mtx", np.array([[0, 0, 0, 1]]), fmt="%.1f") 


            print("###### do some demographic for the paper ###########")

            print(PAZ_data.shape)
            numeric_cols = PAZ_data.select_dtypes(include='number').columns

            print(f"MCI at BL {PAZ_data.groupby(['MCI_T0']).size()}")
            print(f"MCI at BL by group {PAZ_data.groupby('Fenotipo')['MCI_T0'].value_counts()}")
            table_mci = pd.crosstab(PAZ_data['Fenotipo'], PAZ_data['MCI_T0'])
            print("table mci")
            print(table_mci)
            print(f"fisher exact test {fisher_exact(table_mci)}")

            print(f"Sex {PAZ_data.groupby(['Sex']).size()}")
            print(f"Sex by group {PAZ_data.groupby('Fenotipo')['Sex'].value_counts()}")
            table_sex = pd.crosstab(PAZ_data['Fenotipo'], PAZ_data['Sex'])
            print(f"fisher exact test {fisher_exact(table_sex)}")

            print("#####")
            PAZ_data['T1_T0_Diff_Days'] = PAZ_data['T1_T0_Diff_Days'] /365
            for g in ['Eta','Durata malattia (anni)','T1_T0_Diff_Days','UMSARS I','UMSARS I.1']: #,'UMSARS II','UMSARS II.1']:
                print(f"#### {g} ####")
                # all statistics 
                #all_shapiro = shapiro(PAZ_data[g])
                #if all_shapiro.pvalue <= 0.05:
                #    print(colored(f"{g} all group: {PAZ_data[g].median()} {PAZ_data[g].agg(Q1=lambda x: x.quantile(0.25),    Q3=lambda x: x.quantile(0.75), IQR=lambda x: x.quantile(0.75) - x.quantile(0.25))}","red"))
                #else:
                #    print(colored(f"{g} all group: {PAZ_data[g].mean()} {PAZ_data[g].std()}","green"))
                    
                
                # statiscs by group 
                group_P = PAZ_data[PAZ_data['Fenotipo'] == 'P'][g] # PAZ_data[PAZ_data['MCI_T0'] == 0][g] 
                group_C = PAZ_data[PAZ_data['Fenotipo'] == 'C'][g] # PAZ_data[PAZ_data['MCI_T0'] == 1][g]
                p_shapiro = shapiro(group_P)
                c_shapiro = shapiro(group_C)
               
                if p_shapiro.pvalue >= 0.05 and c_shapiro.pvalue >= 0.05:  # if normally distributed
                    print(colored(f"{g} MSAP and MSAC are normall","cyan"))
                    print(colored(f" ttest {ttest_ind(group_P, group_C)}","cyan"))
                    print(colored(f"{g} MEAN by group: {PAZ_data.groupby('Fenotipo')[g].mean()}","cyan"))   
                    print(colored(f"{g} STD by group: {PAZ_data.groupby('Fenotipo')[g].std()}","cyan"))   
                else:
                    print(colored(f"{g} MSAP or MSAC are not normal","magenta"))
                    print(colored(f"MSAP pvalue {p_shapiro.pvalue}","magenta"))
                    print(colored(f"MSAC pvalue {c_shapiro.pvalue}","magenta"))
                    print(colored(f"{g} MEDIAN by group: {PAZ_data.groupby('Fenotipo')[g].median()}","magenta"))
                    print(colored(f"{g} IQR by group: {PAZ_data.groupby('Fenotipo')[g].agg(Q1=lambda x: x.quantile(0.25),    Q3=lambda x: x.quantile(0.75), IQR=lambda x: x.quantile(0.75) - x.quantile(0.25))}","magenta"))   
                    print(colored(f"{g} MEAN by group: {PAZ_data.groupby('Fenotipo')[g].mean()}","magenta"))   
                    print(colored(f"{g} STD by group: {PAZ_data.groupby('Fenotipo')[g].std()}","magenta"))   
                    print(colored(f" mannwithney {mannwhitneyu(group_P, group_C, alternative='two-sided')}","magenta"))
                 



    except Exception as e:
        logger.error(f"An error occurred during the script: {e}")
    finally:
        logger.info("Finished processing data.")

if __name__ == "__main__":
    sys.exit(main(sys.argv))