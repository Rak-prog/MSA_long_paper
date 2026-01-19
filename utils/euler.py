import numpy as np
import pandas as pd 


fsqc_df_sesT0 = pd.read_csv("/mnt/NAS_Progetti/PSP-MSA/daPellecchia/Bids_Salerno/derivatives/freesurfer_qualitycheck_allsubjects/ses-T0/fsqc-results.csv")
fsqc_df_sesT1 = pd.read_csv("/mnt/NAS_Progetti/PSP-MSA/daPellecchia/Bids_Salerno/derivatives/freesurfer_qualitycheck_allsubjects/ses-T1/fsqc-results.csv")
print(fsqc_df_sesT1.shape)
subj_study_long =  ['sub-05_ses-T0', 'sub-07_ses-T0', 'sub-08_ses-T0','sub-09_ses-T0','sub-11_ses-T0',\
                  'sub-12_ses-T0', 'sub-13_ses-T0', 'sub-15_ses-T0','sub-17_ses-T0','sub-19_ses-T0',\
                  'sub-20_ses-T0', 'sub-21_ses-T0', 'sub-23_ses-T0','sub-25_ses-T0','sub-26_ses-T0',\
                  'sub-27_ses-T0', 'sub-28_ses-T0', 'sub-29_ses-T0','sub-30_ses-T0','sub-31_ses-T0',\
                  'sub-34_ses-T0', 'sub-35_ses-T0', 'sub-36_ses-T0','sub-37_ses-T0','sub-39_ses-T0',\
                  'sub-43_ses-T0', 'sub-47_ses-T0', 'sub-48_ses-T0','sub-50_ses-T0','sub-52_ses-T0',\
                  'sub-54_ses-T0', 'sub-55_ses-T0', 'sub-56_ses-T0','sub-57_ses-T0','sub-58_ses-T0','sub-61_ses-T0']


fsqc_df_sesT0['topo_proxy'] = (fsqc_df_sesT0['holes_lh'] + fsqc_df_sesT0['holes_rh']) 
fsqc_df_sesT0['euler'] = np.abs(2 - 2*((fsqc_df_sesT0['defects_lh'] + fsqc_df_sesT0['defects_rh'])))
fsqc_df_sesT0['euler_orig'] = 2 - 2*((fsqc_df_sesT0['defects_lh'] + fsqc_df_sesT0['defects_rh']))
fsqc_df_sesT0['euler_topo'] = 2 - 2*((fsqc_df_sesT0['topo_lh'] + fsqc_df_sesT0['topo_rh']))
fsqc_df_sesT0['euler_holes'] = 2 - 2*((fsqc_df_sesT0['holes_lh'] + fsqc_df_sesT0['holes_rh']))
cols = ['euler_orig', 'euler_topo', 'euler_holes']
corr_matrix = fsqc_df_sesT0[cols].corr(method='pearson')
print(corr_matrix)


fsqc_df_sesT1['euler_orig'] = 2 - 2*((fsqc_df_sesT1['defects_lh'] + fsqc_df_sesT1['defects_rh']))

# drop everything but MSA ses-T0fsqc_df_sesT0
fsqc_df_sesT0 = fsqc_df_sesT0[fsqc_df_sesT0['subject'].str.contains('|'.join(subj_study_long), na=False)]
# to check controls at t0 only 
#fsqc_df_sesT0 = fsqc_df_sesT0[fsqc_df_sesT0['subject'].str.contains('HC', na=False)]

#fsqc_df_sesT1 = fsqc_df_sesT1[fsqc_df_sesT1['subject'].str.contains('|'.join(subj_study_long), na=False)]

print(fsqc_df_sesT0['euler_orig'].sort_values())
print("##############")
print(fsqc_df_sesT1['euler_orig'].sort_values())







