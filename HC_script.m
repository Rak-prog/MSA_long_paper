clear
clc
close all
input_path='/home/riccardo/codici_progetti/Salerno/stats_results/';

%%%% analysis CROSS (with also left and right) between MSA and HC 

%%%%% do_CT and do_VOL are mutually exclusive
do_CT = 0;
do_VOL = 1;
two_way = false;
match_age_sex = false;
atlas_name = 'DKT'; % DKT or Scfdhaefer
PARC = 400;
NNET = 17;
left_right_merge = true; 

alpha = 0.05;
alpha_bonf = 0.0016;
cross_data = 'crossMSA';  % to select all MSA cross sample or long only 
%%%%
if (do_CT + do_VOL) ~= 1
    error('do_CT and do_VOl are mutually exclusive');
end

if do_CT == 1 & do_VOL == 0

    varname = 'Thickness';
    if strcmp(atlas_name,'DKT')
        left  = readtable([input_path, 'thickness_', cross_data, '_lh.txt']);  % usare qui i risultati del cross non del long 
        right = readtable([input_path, 'thickness_', cross_data, '_rh.txt']);
    
        left_HC = readtable([input_path, 'thickness_crossHC_lh.txt']);
        right_HC = readtable([input_path, 'thickness_crossHC_rh.txt']);
    
        right_HC.temporalpole = [];
        if strcmp(cross_data, 'crossMSAALL')
            right.temporalpole = [];
        end 

        % remove patients that are exluded in the long analysis plue sub-59
        % which has no defined phenotype 
        rowsToRemove = contains(left.lh_aparc_DKTatlas_thickness, {'sub-18_ses-T0','sub-45_ses-T0','sub-53_ses-T0','sub-59_ses-T0'});
        left(rowsToRemove,:)=[];
        right(rowsToRemove,:)=[];
        
        % rename columns
        left = renamevars(left,["lh_aparc_DKTatlas_thickness"],["Subject"]);
        right = renamevars(right,["rh_aparc_DKTatlas_thickness"],["Subject"]);

        left_HC = renamevars(left_HC,["lh_aparc_DKTatlas_thickness"],["Subject"]);
        right_HC = renamevars(right_HC,["rh_aparc_DKTatlas_thickness"],["Subject"]);
    
        % check 
        if ~isequal(left.Properties.VariableNames, right.Properties.VariableNames)
        error('Left and right CT variable names do not match.');
        end
        if ~isequal(left_HC.Properties.VariableNames, right_HC.Properties.VariableNames)
            error('Left and right HC CT variable names do not match.');
        end
    elseif strcmp(atlas_name,'Schaefer')
        left  = readtable([input_path, 'thickness_crossMSA_ses-T0_Schaefer_Parc',num2str(PARC),...
                            '_Net',num2str(NNET),'_lh.txt']);  % usare qui i risultati del cross non del long 
        right = readtable([input_path, 'thickness_crossMSA_ses-T0_Schaefer_Parc',num2str(PARC),...
                            '_Net',num2str(NNET),'_rh.txt']);

        left_HC = readtable([input_path, 'thickness_crossHC_ses-T0_Schaefer_Parc',num2str(PARC),...
                            '_Net',num2str(NNET),'_lh.txt']);
        right_HC = readtable([input_path, 'thickness_crossHC_ses-T0_Schaefer_Parc',num2str(PARC),...
                            '_Net',num2str(NNET),'_rh.txt']);
    end 

    left = left(:,1:end-3);  % last three MeanThickness, BrainSegNoVent, eTIv, first one is Subject col
    right = right(:,1:end-3);
    
    left_HC = left_HC(:,1:end-3);
    right_HC = right_HC(:,1:end-3);

    labels=erase(left.Properties.VariableNames,'');

elseif do_CT == 0 & do_VOL == 1
    varname = 'Volume Normalized';

    vol_synthseg = readtable('/home/riccardo/codici_progetti/Salerno/synthseg_results_cross_MSA.csv'); % important load here the longitudinal results of synthseg
    vol_synthseg_HC = readtable('/home/riccardo/codici_progetti/Salerno/synthseg_results_cross_HC.csv'); % load here the cross sectional
    to_drop = {'leftCerebralCortex','rightCerebralCortex','csf','leftCerebralWhiteMatter' 'rightCerebralWhiteMatter',...
        'leftInferiorLateralVentricle','rightInferiorLateralVentricle','leftLateralVentricle','rightLateralVentricle',...
        'rightHippocampus','leftHippocampus','rightAmygdala','leftAmygdala','rightAccumbensArea','leftAccumbensArea',...
        'leftVentralDC','rightVentralDC'};
    vol_synthseg(:,to_drop) = [];
    vol_synthseg_HC(:,to_drop) = [];
    if strcmp(cross_data , 'crossMSAALL') % remove only the outliers and usbjects not used for a total of 66
        rowsToRemove = contains(vol_synthseg.Subject, {'sub-18_ses-T0','sub-45_ses-T0','sub-53_ses-T0','sub-59_ses-T0'}); 
    elseif strcmp(cross_data , 'crossMSA') % retain only the long subjects for the volume only
        subjectsToKeep = {'sub-05_ses-T0', 'sub-07_ses-T0', 'sub-08_ses-T0','sub-09_ses-T0','sub-11_ses-T0',...
                  'sub-12_ses-T0', 'sub-13_ses-T0', 'sub-15_ses-T0','sub-17_ses-T0','sub-19_ses-T0',...
                  'sub-20_ses-T0', 'sub-21_ses-T0', 'sub-23_ses-T0','sub-25_ses-T0','sub-26_ses-T0',...
                  'sub-27_ses-T0', 'sub-28_ses-T0', 'sub-29_ses-T0','sub-30_ses-T0','sub-31_ses-T0', ...
                  'sub-34_ses-T0', 'sub-35_ses-T0', 'sub-36_ses-T0','sub-37_ses-T0','sub-39_ses-T0',...
                  'sub-43_ses-T0', 'sub-47_ses-T0', 'sub-48_ses-T0','sub-50_ses-T0','sub-52_ses-T0',...
                  'sub-54_ses-T0', 'sub-55_ses-T0', 'sub-56_ses-T0','sub-57_ses-T0','sub-58_ses-T0','sub-61_ses-T0'};
        rowsToRemove = ~ismember(vol_synthseg.Subject, subjectsToKeep);
    else
        error('input name is wrong')
    end
    vol_synthseg(rowsToRemove,:)=[];
    % normalize the volumetrics values
    vol_synthseg(:,2:end-1) = array2table(table2array(vol_synthseg(:,2:end-1)) ./ table2array(vol_synthseg(:,1)).*100);
    vol_synthseg_HC(:,2:end-1) = array2table(table2array(vol_synthseg_HC(:,2:end-1)) ./ table2array(vol_synthseg_HC(:,1)).*100);

    left = vol_synthseg(: , contains(vol_synthseg.Properties.VariableNames, 'left') | contains(vol_synthseg.Properties.VariableDescriptions,'Subject') | contains(vol_synthseg.Properties.VariableDescriptions,'brain-stem'));
    right = vol_synthseg(: , contains(vol_synthseg.Properties.VariableNames, 'right') | contains(vol_synthseg.Properties.VariableDescriptions,'Subject') | contains(vol_synthseg.Properties.VariableDescriptions,'brain-stem'));

    left_HC = vol_synthseg_HC(: , contains(vol_synthseg_HC.Properties.VariableNames, 'left') | contains(vol_synthseg_HC.Properties.VariableDescriptions,'Subject') | contains(vol_synthseg_HC.Properties.VariableDescriptions,'brain-stem'));
    right_HC = vol_synthseg_HC(: , contains(vol_synthseg_HC.Properties.VariableNames, 'right') | contains(vol_synthseg_HC.Properties.VariableDescriptions,'Subject') | contains(vol_synthseg_HC.Properties.VariableDescriptions,'brain-stem'));

    left.CerebellumTot = left.leftCerebellumWhiteMatter + left.leftCerebellumCortex;
    right.CerebellumTot = right.rightCerebellumWhiteMatter + right.rightCerebellumCortex;

    left_HC.CerebellumTot = left_HC.leftCerebellumWhiteMatter + left_HC.leftCerebellumCortex;
    right_HC.CerebellumTot = right_HC.rightCerebellumWhiteMatter + right_HC.rightCerebellumCortex;

    left = left(:,[1:7 9 8]); %, left.Properties.VariableNames(end), left.Properties.VariableNames(end-1)];
    right = right(:,[1:7 9 8]);

    left_HC = left_HC(:,[1:7 9 8]); %, left.Properties.VariableNames(end), left.Properties.VariableNames(end-1)];
    right_HC = right_HC(:,[1:7 9 8]);

    left.Properties.VariableNames = erase(left.Properties.VariableNames,'left');
    right.Properties.VariableNames = erase(right.Properties.VariableNames,'right');

    left_HC.Properties.VariableNames = erase(left_HC.Properties.VariableNames,'left');
    right_HC.Properties.VariableNames = erase(right_HC.Properties.VariableNames,'right');

    if isequal(left.Properties.VariableNames , left_HC.Properties.VariableNames) ~= 1; error('this labels must match'); end
    if isequal(right.Properties.VariableNames , right_HC.Properties.VariableNames)~= 1; error('this labels must match'); end
    newOrder = left.Properties.VariableNames;
    right = right(:, newOrder);  % to move the columns to amtch left, right patients and controls
    right_HC = right_HC(:, newOrder);

    if isequal(left.Properties.VariableNames , right.Properties.VariableNames) ~= 1; error('this labels must match'); end
    if isequal(left_HC.Properties.VariableNames , right_HC.Properties.VariableNames) ~= 1; error('this labels must match'); end
    if isequal(right.Properties.VariableNames , right_HC.Properties.VariableNames) ~= 1; error('this labels must match'); end

    vol_brainstem = readtable('/home/riccardo/codici_progetti/Salerno/brainstem_results_cross_MSA.csv');
    vol_brainstem_HC = readtable('/home/riccardo/codici_progetti/Salerno/brainstem_results_cross_HC.csv');

    % normalize both volumes of synthseg and brainstem segmentations
    % join with the ICV calculated with synthseg in order to norm by ICV
    % important: when doing this merge I am loosing already the outlier
    % subject which are removed from the vol_synthseg table :)
    vol_brainstem = innerjoin(vol_brainstem, vol_synthseg(:,{'Subject','totalIntracranial'}), 'Keys',{'Subject'});
    vol_brainstem(:,2:end) = array2table(table2array(vol_brainstem(:,2:end)) ./ table2array(vol_brainstem(:,end)).*100);

    vol_brainstem_HC = innerjoin(vol_brainstem_HC, vol_synthseg_HC(:,{'Subject','totalIntracranial'}), 'Keys',{'Subject'});
    vol_brainstem_HC(:,2:end) = array2table(table2array(vol_brainstem_HC(:,2:end)) ./ table2array(vol_brainstem_HC(:,end)).*100);
end
   
if do_CT == 0 & do_VOL == 1
    left = left(contains(left.Subject, 'ses-T0'),1:end);
    right = right(contains(right.Subject, 'ses-T0'),1:end);
    labels=left.Properties.VariableNames;

    left_HC = left_HC(:,1:end);
    right_HC = right_HC(:,1:end);


    brainstem_T0 = vol_brainstem(contains(vol_brainstem.Subject, 'ses-T0'),1:end-1); % session T0 excluding totalIntracranial
    labels_bs =  vol_brainstem.Properties.VariableNames(2:end); % excluding subject
    brainstem_HC =  vol_brainstem_HC(:,1:end-1);

    % add Pons to volumes, it is unilateral so I add it only to the left 
    if isequal(brainstem_HC.Subject, left_HC.Subject) ~= 1; error('this labels must match'); 
    else 
        left_HC.Pons = brainstem_HC.Pons; 
        left_HC = movevars(left_HC, 'Pons', 'Before', 'brain_stem');
    end
    
    if isequal(brainstem_T0.Subject, left.Subject) ~= 1; error('this labels must match'); 
    else 
        left.Pons = brainstem_T0.Pons; 
        left = movevars(left, 'Pons', 'Before', 'brain_stem');
    end
end

% demographic information
demo = readtable('/home/riccardo/codici_progetti/Salerno/demographic_plus_dicom_processed_final.csv');
if strcmp(unique(demo.Session), 'ses-T0')
    demo.Subj_ID = strcat(demo.Subj_ID ,'_', demo.Session);
    demo = renamevars(demo,["Subj_ID"],["Subject"]);
end

% merge left and right
if do_CT == 1 & do_VOL == 0
    if strcmp(atlas_name, 'DKT')
        areas = {'superiorfrontal', 'caudalmiddlefrontal', 'rostralmiddlefrontal','parstriangularis','parsopercularis','parsorbitalis','lateralorbitofrontal','medialorbitofrontal','precentral','paracentral',...
                'superiorparietal', 'inferiorparietal', 'supramarginal','postcentral','precuneus',...
                'superiortemporal','middletemporal','inferiortemporal', 'fusiform', 'transversetemporal','parahippocampal','entorhinal',...
                'lateraloccipital','lingual','cuneus','pericalcarine',...
                'rostralanteriorcingulate', 'caudalanteriorcingulate',...
                'posteriorcingulate','isthmuscingulate',...
                'insula'};
        tot_T0 = (left(:,2:end) + right(:,2:end))./2;
        tot_HC = (left_HC(:,2:end) + right_HC(:,2:end))./2;
        tot_T0.Subject = left.Subject;
        tot_HC.Subject = left_HC.Subject; 
    elseif strcmp(atlas_name , 'Schaefer')
         if NNET == 7
            nets_name = ["Vis","SomMot","DorsAttn","SalVentAttn","Limbic","Cont","Default"];
            if left_right_merge == false
                areas = {'LH_Vis','LH_SomMot','LH_DorsAttn','LH_SalVentAttn','LH_Limbic','LH_Cont','LH_Default',...
                        'RH_Vis','RH_SomMot','RH_DorsAttn','RH_SalVentAttn','RH_Limbic','RH_Cont','RH_Default'};
            elseif left_right_merge == true
                areas = {'Vis','SomMot','DorsAttn','SalVentAttn','Limbic','Cont','Default'};
            end
        elseif NNET == 17
            nets_name = ["VisCent","VisPeri","SomMotA","SomMotB","DorsAttnA","DorsAttnB","SalVentAttnA",...
                         "SalVentAttnB","LimbicA","LimbicB","ContA","ContB","ContC",...
                         "DefaultA","DefaultB","DefaultC","TempPar"];
            if left_right_merge == false
                areas = {'LH_VisCent','LH_VisPeri','LH_SomMotA','LH_SomMotB','LH_DorsAttnA','LH_DorsAttnB','LH_SalVentAttnA',...
                     'LH_SalVentAttnB','LH_LimbicA','LH_LimbicB','LH_ContA','LH_ContB','LH_ContC',...
                     'LH_DefaultA','LH_DefaultB','LH_DefaultC','LH_TempPar',...
                     'RH_VisCent','RH_VisPeri','RH_SomMotA','RH_SomMotB','RH_DorsAttnA','RH_DorsAttnB','RH_SalVentAttnA',...
                     'RH_SalVentAttnB','RH_LimbicA','RH_LimbicB','RH_ContA','RH_ContB','RH_ContC',...
                     'RH_DefaultA','RH_DefaultB','RH_DefaultC','RH_TempPar'};
            elseif left_right_merge == true
                areas = {'VisCent','VisPeri','SomMotA','SomMotB','DorsAttnA','DorsAttnB','SalVentAttnA',...
                     'SalVentAttnB','LimbicA','LimbicB','ContA','ContB','ContC',...
                     'DefaultA','DefaultB','DefaultC','TempPar'};
            end
         end
        left.Properties.VariableNames = erase(left.Properties.VariableNames, 'lh_');
        right.Properties.VariableNames = erase(right.Properties.VariableNames, 'rh_');
        left_HC.Properties.VariableNames = erase(left_HC.Properties.VariableNames, 'lh_');
        right_HC.Properties.VariableNames = erase(right_HC.Properties.VariableNames, 'rh_');

        keyname = ['Schaefer2018_', num2str(PARC),'Parcels_', num2str(NNET),'Networks_order_thickness'];
   
        if isequal(right.(keyname) , left.(keyname))
             data = innerjoin(left, right, 'Keys' , keyname);
        else
            error('Subject id are not matching, please double check')
        end
        if isequal(right_HC.(keyname) , left_HC.(keyname))
             data_HC = innerjoin(left_HC, right_HC, 'Keys' , keyname);
        else
            error('Subject id are not matching, please double check')
        end 
        tot_T0 = merge_networks(data(:,2:end), nets_name, left_right_merge);
        tot_HC = merge_networks(data_HC(:,2:end),nets_name, left_right_merge);
        tot_T0.Subject = data.(keyname);
        tot_HC.Subject = data_HC.(keyname);
        left.Subject = left.(keyname);
        right.Subject = right.(keyname);
        left_HC.Subject = left_HC.(keyname);
        right_HC.Subject = right_HC.(keyname);
    end
elseif do_CT == 0 & do_VOL == 1
     areas = {'Putamen', 'CerebellumWhiteMatter','CerebellumCortex','Pons'};
     pons = left.Pons;
     pons_HC = left_HC.Pons; 

     left_HC.Pons = []; % togli il pons per il calcolo, poi lo rimetti 
     left.Pons = [];
     tot_T0  = (left(:,1:end-1) + right(:,1:end-1));
     tot_HC  = (left_HC(:,1:end-1) + right_HC(:,1:end-1));
     
     left.Pons = pons; % re-insert pons
     left_HC.Pons = pons_HC;
     
     tot_T0.Pons = pons;
     tot_HC.Pons = pons_HC;
     
     tot_T0.Subject = left.Subject;
     tot_HC.Subject = left_HC.Subject;
end

% do an anova correcting for age and sex, so merge demo with the data 
left_demo = innerjoin( [left ; left_HC ], demo, 'Keys',{'Subject'});
right_demo = innerjoin( [right ; right_HC ], demo, 'Keys',{'Subject'});
tot_demo = innerjoin( [tot_T0 ; tot_HC ], demo, 'Keys',{'Subject'});
cols_name = {'Subject','GROUP','sex_final'};

for col = 1:numel(cols_name)
    left_demo.(cols_name{col}) = categorical(left_demo.(cols_name{col}));
    right_demo.(cols_name{col}) = categorical(right_demo.(cols_name{col}));
    tot_demo.(cols_name{col}) = categorical(tot_demo.(cols_name{col}));
end

% 3 way ANOVA in enabled
if two_way == true
    phenotype_modified = tot_demo.Phenotype;
    tot_demo.GROUP(strcmp(tot_demo.Phenotype , 'P')) = 'MSAP';
    tot_demo.GROUP(strcmp(tot_demo.Phenotype , 'C')) = 'MSAC';
    tot_demo.GROUP = categorical( tot_demo.GROUP);
    tot_demo.GROUP = removecats(tot_demo.GROUP, 'MSA');
    % Create the GROUP variable
end

if match_age_sex == false
    [recap_all_2groups,relative_all, realtive_msap, realtive_msac, CI_high, partial_eta] = do_between_group_anova(tot_demo,areas,alpha, false, two_way);
    
   
    % adding partial eta (effect size to the table)
    for i = 1:numel(areas)
        roi = areas{i};
        eta_colname = [roi '_eta'];
        pcol_idx = find(strcmp(recap_all_2groups.Properties.VariableNames, roi));
    
        if isempty(pcol_idx)
            error('Column %s not found in recap_table', roi);
        end
    
        eta_col = repmat(partial_eta(i), height(recap_all_2groups), 1);
    
        % Insert eta column immediately AFTER the p-value column
        recap_all_2groups = addvars(recap_all_2groups, ...
            eta_col, ...
            'After', pcol_idx, ...
            'NewVariableNames', eta_colname);
    end
    if two_way == false % save only the single group results (no MSA p ands MSA C group effect)
        if do_CT == 1 && do_VOL == 0
           if strcmp(atlas_name , 'DKT') == 1
               writetable(tot_demo, ['/home/riccardo/codici_progetti/MSA_paper/MSA_long_paper/matrices/CT_cross_HC_vs_', cross_data , '.csv']);
               writetable(recap_all_2groups,[ './matrices/CT_recap_HC_vs_',cross_data,'.tsv'], 'FileType', 'text', 'Delimiter', '\t');
               writetable(recap_all_2groups,[ './matrices/CT_recap_HC_vs_',cross_data,'.csv']);
           elseif strcmp(atlas_name, 'Schaefer') == 1
               writetable(tot_demo, ['/home/riccardo/codici_progetti/MSA_paper/MSA_long_paper/matrices/CT_cross_HC_vs',cross_data,'_Schaefer_Parc',num2str(PARC),'_NET', num2str(NNET),'.csv']);
               writetable(recap_all_2groups,[ '/home/riccardo/codici_progetti/MSA_paper/MSA_long_paper/matrices/CT_recap_HC_vs',cross_data,'_Schaefer_Parc',num2str(PARC),'_NET', num2str(NNET),'.tsv'], 'FileType', 'text', 'Delimiter', '\t');
               writetable(recap_all_2groups,[ '/home/riccardo/codici_progetti/MSA_paper/MSA_long_paper/matrices/CT_recap_HC_vs',cross_data,'_Schaefer_Parc',num2str(PARC),'_NET', num2str(NNET),'.csv']);
           end
        elseif do_CT == 0 && do_VOL == 1
           writetable(tot_demo, ['/home/riccardo/codici_progetti/MSA_paper/MSA_long_paper/matrices/VOL_cross_HC_vs_', cross_data , '.csv']);
           writetable(recap_all_2groups,[ './matrices/VOL_recap_HC_vs_',cross_data,'.tsv'], 'FileType', 'text', 'Delimiter', '\t');
           writetable(recap_all_2groups,[ './matrices/VOL_recap_HC_vs_',cross_data,'.csv']);
        end 
    end
    if two_way == true
    
        row_group_msap = strcmp(recap_all_2groups.Labels, 'GROUP_MSAP');
        row_group_msac = strcmp(recap_all_2groups.Labels, 'GROUP_MSAC');
        pval_family_msap = table2array(recap_all_2groups(row_group_msap,3:3:end)); 
        pval_family_msac = table2array(recap_all_2groups(row_group_msac,3:3:end)); 
        hh_msap = fdr_bh(pval_family_msap,alpha,'pdep','yes');
        hh_msac = fdr_bh(pval_family_msac,alpha,'pdep','yes');
    elseif two_way == false
    
        row_group = strcmp(recap_all_2groups.Labels, 'GROUP_MSA');
        pval_family = table2array(recap_all_2groups(row_group,3:3:end)); 
        beta_family = recap_all_2groups(row_group,2:3:end); 
        pval_family_table = (recap_all_2groups(row_group,3:3:end));
        hh = fdr_bh(pval_family,alpha,'pdep','yes');
    
        if sum(hh)>=1
            significant = pval_family_table(: , hh)
            beta_family(:,hh)
            CI_low(row_group, hh)
            CI_high(row_group, hh)
            partial_eta(hh)
        else
            disp('No significant ROI after Benjiamini-Hochberg correction')
        end
    
    end 
    if two_way == false % save only the single group results 
        if do_CT == 1
            if strcmp(atlas_name , 'DKT') == 1
                writetable(recap_all_2groups, ['/home/riccardo/codici_progetti/MSA_paper/MSA_long_paper/matrices/anova_tot_ct_cross_2groups.csv']);
            elseif strcmp(atlas_name , 'Schaefer') == 1
                writetable(recap_all_2groups, ['/home/riccardo/codici_progetti/MSA_paper/MSA_long_paper/matrices/anova_tot_ct_cross_2groups','_Schaefer_Parc',num2str(PARC),'_NET', num2str(NNET),'.csv']);
            end
        elseif do_VOL == 1
            writetable(recap_all_2groups, ['/home/riccardo/codici_progetti/MSA_paper/MSA_long_paper/matrices/anova_tot_vol_cross_2groups.csv']);
        end
    end
elseif match_age_sex == true

    disp('doing sex and age matched analysis')
    % filter tot_demo to get only HC age and sex matched
    MSA = tot_demo(ismember(tot_demo.GROUP, {'MSA'}), :);
    HC  = tot_demo(ismember(tot_demo.GROUP, {'HC'}),  :);

    nMSA = height(MSA);
    
    matched_HC = table();     % output
    used_HC = false(height(HC),1);  % flag to avoid reuse
    
    age_tol = 1; 


    for i = 1:nMSA

        age_i = MSA.age_final(i);
        sex_i = MSA.sex_final(i);

        % candidate HC: same sex, age within tolerance, not used
        idx = ~used_HC & ...
            HC.sex_final == sex_i & ...
            abs(HC.age_final - age_i) <= age_tol;

        candidates = HC(idx,:);

        if isempty(candidates)
            warning('No HC match for MSA subject %d', i);
            continue
        end

        % choose closest in age
        [~, k] = min(abs(candidates.age_final - age_i));
        selected = candidates(k,:);

        % store
        matched_HC = [matched_HC; selected];

        % mark as used
        hc_idx_global = find(idx);
        used_HC(hc_idx_global(k)) = true;
    end

    tot_matched = [MSA; matched_HC];

    % check 
    sex_hc = table2array(tot_matched(ismember(tot_matched.GROUP, {'HC'}), 'sex_final'));
    sex_msa = table2array(tot_matched(ismember(tot_matched.GROUP, {'MSA'}), 'sex_final'));
    disp(countcats(sex_hc))
    disp(countcats(sex_msa))

    age_hc = table2array(tot_matched(ismember(tot_matched.GROUP, {'HC'}), 'age_final'));
    age_msa = table2array(tot_matched(ismember(tot_matched.GROUP, {'MSA'}), 'age_final'));
    disp(['avg and sd of hc ', num2str(mean(age_hc)), '  ' , num2str(std(age_hc))])
    disp(['avg and sd of MSA ', num2str(mean(age_msa)),'  ', num2str(std(age_msa))])

    h_diff = NaN*ones(numel(areas),1);
    p_diff = NaN*ones(numel(areas),1);
    for roi = 1:numel(areas)
        roi_hc = table2array(tot_matched(ismember(tot_matched.GROUP, {'HC'}), areas(roi)));
        roi_msa = table2array(tot_matched(ismember(tot_matched.GROUP, {'MSA'}), areas(roi)));

        % --- Normality (per group) ---
        hx = lillietest(roi_hc, 'Alpha', alpha); % 0 => normal-like
        hy = lillietest(roi_msa, 'Alpha', alpha);
        
        is_normal = (hx == 0) && (hy == 0);
        if is_normal
            % Welch t-test (unequal variances)
            [h_diff(roi),p_diff(roi),ci,stats] = ttest2(roi_hc, roi_msa, 'Vartype','unequal', 'Alpha', alpha);
            out.test = "Welch t-test (ttest2, Vartype=unequal)";
            out.statistic = stats.tstat;
            out.df = stats.df;
        else
            [p_diff(roi),h_diff(roi),stats] = ranksum(roi_hc, roi_msa, 'Alpha', alpha);
            out.test = "Mann–Whitney (ranksum)";
            out.statistic = stats.ranksum;
            out.df = NaN;
            ci = [NaN NaN];
        end
    end
    % FDR 
    fdr_bh(p_diff,alpha,'pdep','yes')
end


%%% check with a plot effect sizes 

vec = 1:numel(areas);

figure; hold on; box off;
plot(vec, effect_size_2groups,'o')
plot(vec(hh==1), effect_size_2groups(hh==1),'r.')
plot(vec , partial_eta, 'g*')
legend('f square','f sqaure in FDR corrected ROI','partial eta')
set(gca,'XTick',vec,'XTickLabel',areas)

%%% correlation quickly on the cross
for rr = 1:numel(areas)
    
    lm_dd = fitlm(tot_demo, ['Disease_Duration', ' ~  age_final + sex_final + ', areas{rr}]);
    lm_ums = fitlm(tot_demo, ['UMSARSI', ' ~ age_final + sex_final + ', areas{rr}]);
    
    if lm_dd.Coefficients.pValue(2) <= alpha
        disp(['Disease duration significant for ', areas{rr}, ' at alpha ', num2str(alpha), ...
             ' with p value ', num2str(round(lm_dd.Coefficients.pValue(2),5)), ...
             ' with beta estimate ', num2str(round(lm_dd.Coefficients.Estimate(2),5))]);
    end

    if lm_ums.Coefficients.pValue(2) <= alpha
        disp(['UMSARS significant for ', areas{rr}, ' at alpha ', num2str(alpha), ...
             ' with p value ', num2str(round(lm_ums.Coefficients.pValue(2),5)), ...
             ' with beta estimate ', num2str(round(lm_ums.Coefficients.Estimate(2),5))]);
    end
end

