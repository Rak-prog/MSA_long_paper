clear
clc
close all

load('./matrices/tot_stats_ct_singleparcel.mat')
group = [1
1
2
2
3
3
3
1
1
1
1
3
3
1
1
3
3
3
1
3
2
1
1
1
1
1
3
1
1
3
1
1
1
3
1
3];

age  = table2array(tot_stats(:,33));
sex = table2array(tot_stats(:,34));
n = grp2idx(sex)-1;

data=table2array(tot_stats(:,1:31));

delta = data(1:36,:)-data(37:72,:);

for i = 1 : size(delta,2)
    tbl = table(age(1:36),sex(1:36),group,delta(:,i),'VariableNames',{'age','sex','group','delta'});

    lm = fitlm(tbl,'delta~age+sex+group')
    P(i)=lm.Coefficients{4,4};
end
  
[h, ~, ~]=fdr_bh(P,0.05,'pdep','yes')

subj_to_keep = {'sub-05_ses-T0', 'sub-07_ses-T0', 'sub-08_ses-T0', 'sub-09_ses-T0', 'sub-11_ses-T0', 'sub-12_ses-T0', 'sub-13_ses-T0', ...
                   'sub-15_ses-T0', 'sub-17_ses-T0', 'sub-19_ses-T0', 'sub-20_ses-T0', 'sub-21_ses-T0', 'sub-23_ses-T0', 'sub-25_ses-T0', ...
                   'sub-26_ses-T0', 'sub-27_ses-T0', 'sub-28_ses-T0', 'sub-29_ses-T0', 'sub-30_ses-T0', 'sub-31_ses-T0', 'sub-34_ses-T0', ...
                   'sub-35_ses-T0', 'sub-36_ses-T0', 'sub-37_ses-T0', 'sub-39_ses-T0', 'sub-43_ses-T0', 'sub-47_ses-T0', 'sub-48_ses-T0', ...
                   'sub-50_ses-T0', 'sub-52_ses-T0', 'sub-54_ses-T0', 'sub-55_ses-T0', 'sub-56_ses-T0', 'sub-57_ses-T0', 'sub-58_ses-T0', ...
                   'sub-61_ses-T0'};

[num,txt,raw] = xlsread('./matrices/demographic_plus_dicom_processed_final_analysis_against_combat.xlsx',5);

subj_col = raw(2:end,1);
keep_mask = ismember(subj_col, subj_to_keep);

num = num(keep_mask, :);
group = num(:,1);
sex = num(:,2);
age=num(:,3);
for i = 4:34
    
   tbl = table(age,sex,group,num(:,i),'VariableNames',{'age','sex','group','delta'});
   tbl.sex = categorical(tbl.sex);
   lm = fitlm(tbl,'delta~age+sex+group')
    
    P(i-3)=lm.Coefficients{4,4}; 
    
end
[h, ~, ~]=fdr_bh(P,0.05)

%%
[num,txt,raw] = xlsread('./matrices/demographic_plus_dicom_processed_final_analysis_against_combat.xlsx',6);
subj_col = raw(2:end,1);
keep_mask = ismember(subj_col, subj_to_keep);

num = num(keep_mask, :);
group = num(:,1);
sex = num(:,2);
age=num(:,3);
for i = 4:34
    
   tbl = table(age,sex,group,num(:,i),'VariableNames',{'age','sex','group','delta'});
    tbl.sex = categorical(tbl.sex);
   lm = fitlm(tbl,'delta~age+sex+group')
   P(i-3)=lm.Coefficients{4,4}; 
    
end
[h, ~, ~]=fdr_bh(P,0.05)
