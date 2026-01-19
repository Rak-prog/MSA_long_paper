clear
close all
clc

do_CT = 0;
do_VOL = 1;
atlas_name = 'Schaefer';
NPARC = 400;
NNET = 17;
left_right_merge = true;
fusing_roi=0;
alpha=0.05;
alpha_bonf = 0.0016;
remove_outliers = 0;
diagnostics = 1;
do_lme_correlation = 1;
nsubj = 36; 
 

if do_VOL == 1 & do_CT == 0
    load('./matrices/tot_stats_vol.mat');
    load('./matrices/tot_stats_brainstem_vol.mat');
    if isequal(tot_stats.Subj_ID,tot_brainstem.Subj_ID) & isequal(tot_stats.SottotipoMotorio,tot_brainstem.SottotipoMotorio)
        data = tot_stats;
        data.Pons = tot_brainstem.Pons;
    end
    areas = {'Putamen', 'CerebellumWhiteMatter','CerebellumCortex','Pons'};
    areas_significant = areas;% the areas that survived correction for multiple comparison for the volume at p < 0.001
elseif   do_VOL ==  0 & do_CT == 1
    if strcmp(atlas_name, 'DKT') == 1
        if fusing_roi == 1
            tot_stats = load('./matrices/tot_stats_ct_fused.mat');
            data = tot_stats.tot_stats;
            areas = {'Frontal', 'Parietal', 'Temporal','Occipital','ACC','PCC','Insula'};
        elseif fusing_roi == 0
            tot_stats = load('./matrices/tot_stats_ct_singleparcel.mat');
            data = tot_stats.tot_stats;
            areas = {'superiorfrontal', 'caudalmiddlefrontal', 'rostralmiddlefrontal','parstriangularis','parsopercularis','parsorbitalis','lateralorbitofrontal','medialorbitofrontal','precentral','paracentral',...
                'superiorparietal', 'inferiorparietal', 'supramarginal','postcentral','precuneus',...
                'superiortemporal','middletemporal','inferiortemporal', 'fusiform', 'transversetemporal','parahippocampal','entorhinal',...
                'lateraloccipital','lingual','cuneus','pericalcarine',...
                'rostralanteriorcingulate', 'caudalanteriorcingulate',...
                'posteriorcingulate','isthmuscingulate',...
                'insula'};
           areas_significant = {'superiorfrontal','caudalmiddlefrontal','precentral',...
                              'posteriorcingulate','isthmuscingulate','insula'};
        end
    elseif strcmp(atlas_name, 'Schaefer') == 1
        % process Schaefer results
        atlas_data_left_T0 = readtable(['/home/riccardo/codici_progetti/Salerno/stats_results/thickness_long_ses-T0_Schaefer_Parc', num2str(NPARC),'_Net',num2str(NNET),'_lh.txt']);
        atlas_data_right_T0 = readtable(['/home/riccardo/codici_progetti/Salerno/stats_results/thickness_long_ses-T0_Schaefer_Parc', num2str(NPARC),'_Net',num2str(NNET),'_rh.txt']);
        atlas_data_left_T1 = readtable(['/home/riccardo/codici_progetti/Salerno/stats_results/thickness_long_ses-T1_Schaefer_Parc', num2str(NPARC),'_Net',num2str(NNET),'_lh.txt']);
        atlas_data_right_T1 = readtable(['/home/riccardo/codici_progetti/Salerno/stats_results/thickness_long_ses-T1_Schaefer_Parc', num2str(NPARC),'_Net',num2str(NNET),'_rh.txt']);
        % remove from the tables the column not used
        col2remove = ["Background_FreeSurfer_Defined_Medial_Wall_thickness", "BrainSegVolNotVent", "eTIV"];
        tables = {atlas_data_left_T0, atlas_data_right_T0, atlas_data_left_T1, atlas_data_right_T1};        
        for i = 1:numel(tables)
            t = tables{i};
            t(:, contains(t.Properties.VariableNames, col2remove)) = [];
            tables{i} = t;
        end

        [atlas_data_left_T0, atlas_data_right_T0, atlas_data_left_T1, atlas_data_right_T1] = tables{:};
        atlas_data_left_T0.Properties.VariableNames = erase(atlas_data_left_T0.Properties.VariableNames, 'lh_');
        atlas_data_right_T0.Properties.VariableNames = erase(atlas_data_right_T0.Properties.VariableNames, 'rh_');
        atlas_data_left_T1.Properties.VariableNames = erase(atlas_data_left_T1.Properties.VariableNames, 'lh_');
        atlas_data_right_T1.Properties.VariableNames = erase(atlas_data_right_T1.Properties.VariableNames, 'rh_');
        
        keyname = ['Schaefer2018_', num2str(NPARC),'Parcels_', num2str(NNET),'Networks_order_thickness'];
        data_T0 = innerjoin(atlas_data_left_T0, atlas_data_right_T0, 'Keys' , keyname);
        data_T1 = innerjoin(atlas_data_left_T1, atlas_data_right_T1, 'Keys' , keyname);
        
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
        data_T0_bynet = merge_networks(data_T0(:,2:end), nets_name, left_right_merge);
        data_T1_bynet = merge_networks(data_T1(:,2:end),nets_name, left_right_merge);

        if isequal(data_T0_bynet.Properties.VariableNames , data_T1_bynet.Properties.VariableNames)
            data_T0_bynet.Subj_ID = data_T0.(['Schaefer2018_', num2str(NPARC),'Parcels_', num2str(NNET),'Networks_order_thickness']);
            data_T1_bynet.Subj_ID = data_T1.(['Schaefer2018_', num2str(NPARC),'Parcels_', num2str(NNET),'Networks_order_thickness']);
            data_T0_bynet.Time = repmat("T0", height(data_T0_bynet), 1);
            data_T1_bynet.Time = repmat("T1", height(data_T1_bynet), 1);
            
            id0 = string(data_T0_bynet.Subj_ID);
            id1 = string(data_T1_bynet.Subj_ID);
            
            % remove session + long suffix (adjust if your exact suffix differs)
            id0 = erase(id0, "_ses-T0.long");
            id1 = erase(id1, "_ses-T1.long");
            
            % (optional) also remove plain session if present without .long
            id0 = erase(id0, "_ses-T0");
            id1 = erase(id1, "_ses-T1");

            if isequal(id0 , id1)
                data_schaefer = [data_T0_bynet ;data_T1_bynet];
            end

        else
            error("T0 and T1 by-network tables have different variable names. Check merge_networks output.");
        end

        % merge with demo from the other CT file 
        tot_stats = load('./matrices/tot_stats_ct_singleparcel.mat');
        tot_stats = tot_stats.tot_stats;
        % fix namings 
        data_schaefer.Subj_ID = erase(data_schaefer.Subj_ID, "_ses-T1.long");  
        data_schaefer.Subj_ID = erase(data_schaefer.Subj_ID, "_ses-T0.long"); 
        if isequal(tot_stats.Subj_ID , data_schaefer.Subj_ID)
            data = [data_schaefer , tot_stats(:,{'Et_','Sesso','SottotipoMotorio','MCI_s__1_no_0_',...
                                                  'DurataMalattia_anni_','UMSARSI','UMSARSI_1','UMSARSIcov',...
                                                  'T1_T0_Diff_Days'})];
        else
            error('Subjects IDs are not matching. Please double check')
        end
                
        data.Subj_ID = categorical(data.Subj_ID);
        data.Sesso = categorical(data.Sesso);
        data.Time = categorical(data.Time);
        data.SottotipoMotorio = categorical(data.SottotipoMotorio);
        data.MCI_s__1_no_0_ = categorical(data.MCI_s__1_no_0_);

        areas_significant = {'SomMotA','SalVentAttnA','SalVentAttnB','ContA','ContB','ContC',...
                                'DefaultA','DefaultB','DefaultC'};

    else
        error('Filename is not valid');
    end   
end
% Time * MCI_s__1_no_0_ + Time * SottotipoMotorio + Et_ + Sesso + (1|Subj_ID)
lme_list = {'~ Time * MCI_s__1_no_0_ + Time * SottotipoMotorio + Et_ + Sesso + (1|Subj_ID)'};
lme_list_savename={'total_model'};

if remove_outliers == 1
    nsubj_orig = size(data,1);
    [data_clean , outliers_subj] = insert_nan_where_outlier(data, areas); % removing outliers across time points and group (one unique distribution)
    rowsToRemove = contains(data.Subj_ID, outliers_subj);
    data(rowsToRemove,:)=[];
    if (nsubj_orig - sum(rowsToRemove)) ~= size(data,1)
        error('outliers removal gone wrong')
    end
end

% convert Subj_ID to categorical
data.Subj_ID = categorical(data.Subj_ID);

for type_of_lme = 1:numel(lme_list)
    for roi = 1:numel(areas)
        % loop over type of lme
        areas{roi}

        lme = fitlme(data, [areas{roi},lme_list{type_of_lme}]);
    
        if any(vifs(lme)>=5)
            err(['there is a collinearity issue for model UMSARSI and area ', areas{roi}])
        end

        if roi == 1
            recap_table = table( lme.CoefficientNames(1:end)', 'VariableNames', {'Labels'});
        end

        recap_table.([areas{roi},'_estimate']) = round(lme.Coefficients.Estimate,4);
        recap_table.(areas{roi}) = round(lme.Coefficients.pValue,12);

        % diagnostics of the model if enabled
        if diagnostics == 1
            
            res_raw  = residuals(lme, 'ResidualType', 'Raw');
            res_std  = residuals(lme, 'ResidualType', 'Standardized');
            yhat_fix = fitted(lme, 'Conditional', false);  % fixed-effects fit

            figure;
            subplot(131)
            plotResiduals(lme, 'fitted')
            axis square

            subplot(132)
            qqplot(res_std);
            title(['LME: Q-Q plot of standardized residuals', areas{roi}]);
            axis square
            %plotResiduals(lme, 'probability')

            res_std_centered = (res_std - mean(res_std)) / std(res_std);
            [h_ks, p_ks] = kstest(res_std_centered);
            disp([h_ks, p_ks]);
            if h_ks ~= 0
                disp(['KS TEST IS SIGNFICANT FOR ROI ', areas{roi}])
            end

            F = fitted(lme);
            R = response(lme);
            subplot(133)
            plot(R,F,'rx')
            xlabel('Response')
            ylabel('Fitted')
            title(['LME: Observed Response vs Fitted Response plot ', areas{roi}]);
            axis square
            close all;
        end
    end

   
    disp(['fdr correction for ' , lme_list{type_of_lme}]);
    if strcmp(atlas_name , 'Schaefer') & do_CT == 1
        interaction_over_time_MCI = strcmp(recap_table.Labels, 'Time_T1:MCI_s__1_no_0__1');
        interaction_over_time_Pheno = strcmp(recap_table.Labels, 'Time_T1:SottotipoMotorio_P');
    elseif strcmp(atlas_name , 'DKT') & do_CT == 1
        interaction_over_time_MCI = strcmp(recap_table.Labels, 'MCI_s__1_no_0__1:Time_T1');
        interaction_over_time_Pheno = strcmp(recap_table.Labels, 'SottotipoMotorio_P:Time_T1');
    elseif do_VOL == 1
        interaction_over_time_MCI = strcmp(recap_table.Labels, 'MCI_s__1_no_0__1:Time_T1');
        interaction_over_time_Pheno = strcmp(recap_table.Labels, 'SottotipoMotorio_P:Time_T1');
    end
    [h_lme_interaction_time_MCI, ~, ~] = fdr_bh(table2array(recap_table(interaction_over_time_MCI,3:2:end)),alpha,'pdep','yes')
    [h_lme_interaction_time_Pheno, ~, ~] = fdr_bh(table2array(recap_table(interaction_over_time_Pheno,3:2:end)),alpha,'pdep','yes')

    h_lme_tot_bonf = table2array(recap_table(:,3:2:end)) <= alpha_bonf
    disp('detect if there are more then 1 signficiant values in each roi')
    sum(h_lme_tot_bonf(2:end,:))
    writetable(recap_table,[ 'recap_MSA_long_interaction.tsv'], 'FileType', 'text', 'Delimiter', '\t');
end

% % visualize BL vs FUP
% data.Subj_ID = string(data.Subj_ID); % converting subj_id column back to cell
% dataBL = data(data.Time == 'T0', areas);
% dataFUP = data(data.Time == 'T1', areas);
% subj_listBL = data(data.Time == 'T0', 'Subj_ID');
% subj_listFUP = data(data.Time == 'T1', 'Subj_ID');
% min_val = min([ min(min(table2array(dataBL)))   min(min(table2array(dataFUP)))]);
% max_val = max([ max(max(table2array(dataBL)))   max(max(table2array(dataFUP)))]);
% if do_CT == 0 & do_VOL == 1
%     [hvols,~,~,~,~,increasing_subjects]=visualize(table2array(dataBL), table2array(dataFUP), 'BL vs FUP', [min_val max_val], areas,...
%         'VOL [%]', table2array(subj_listBL), table2array(subj_listFUP), 1, 'b', 'c')
%     exportgraphics(hvols,['/home/riccardo/codici_progetti/Salerno/figure_paper/BLvsFUP_vol_remove_outliers_', num2str(remove_outliers),'.png'], 'Resolution', 300);
% elseif do_CT == 1 & do_VOL == 0
%     [hct,~,~,~,~,increasing_subjects]=visualize(table2array(dataBL), table2array(dataFUP), 'BL vs FUP', [min_val max_val], areas,...
%         'CT [mm]', table2array(subj_listBL), table2array(subj_listFUP), 1, 'b', 'c')
%     exportgraphics(hct,['/home/riccardo/codici_progetti/Salerno/figure_paper/BLvsFUP_ct_fusing_roi_',num2str(fusing_roi),'_remove_outliers_',num2str(remove_outliers),'.png'], 'Resolution', 300);
% end

% % outliers list
% max_len = max(cellfun(@numel, increasing_subjects));
% 
% % Step 2: Pad shorter cells with empty strings
% padded = cell(numel(increasing_subjects), max_len);
% for i = 1:numel(increasing_subjects)
%     row = increasing_subjects{i};
%     padded(i, 1:numel(row)) = row;
% end

% % Step 3: Convert to table
% T = cell2table(padded);
% %T.Properties.VariableNames = repmat('', 1, max_len);
% T.Properties.RowNames = areas;
% 
% % Display the table
% if do_CT == 0 & do_VOL == 1
%     writetable(T,['/home/riccardo/codici_progetti/Salerno/figure_paper/VOL_subj_increase_remove_outliers_',num2str(remove_outliers),'.csv'],'WriteRowNames', true);
% elseif do_CT == 1 & do_VOL == 0
%     writetable(T,['/home/riccardo/codici_progetti/Salerno/figure_paper/CT_subj_increase_fusing_',num2str(fusing_roi),'_remove_outliers_',num2str(remove_outliers),'.csv'],'WriteRowNames', true);
% end

% allCells = table2cell(T);
% % Flatten into a single column
% flatList = allCells(:);

% % Remove empty entries (optional)
% flatList = flatList(~cellfun(@isempty, flatList));
% 
% % Get unique strings
% uniqueStrings = unique(flatList);
% clear lme areas
corr_variables = {'disease_duration_centered','UMSARSIcov_centered'};

% add time at FUP to disease duration at FUP
data.DurataMalattia_anni_adjusted(data.Time=='T1') = data.DurataMalattia_anni_(data.Time=='T0') + data.T1_T0_Diff_Days(data.Time=='T0')./365; 
data.DurataMalattia_anni_adjusted(data.Time=='T0') = data.DurataMalattia_anni_(data.Time=='T0'); 
data.DurataMalattia_anni_adjusted = data.DurataMalattia_anni_adjusted - mean(data.DurataMalattia_anni_adjusted);
if  ~isequal(round(data.DurataMalattia_anni_adjusted(data.Time=='T1') - data.DurataMalattia_anni_adjusted(data.Time=='T0'),4) , round(data.T1_T0_Diff_Days(data.Time=='T0')./365,4))
    error('adjusting disease duration is wrong')
end

data.UMSARSIcov_centered = data.UMSARSIcov - mean(data.UMSARSIcov);

if ~isequal(data.Subj_ID(1:nsubj), data.Subj_ID(nsubj+1:2*nsubj))
    error('Subject order mismatch between T0 and T1');
end

% check that baseline covariates are identical at both timepoints
vars_to_check = {'Et_','Sesso','MCI_s__1_no_0_','SottotipoMotorio'};
for v = vars_to_check
    vname = v{1};
    if any(data.(vname)(1:nsubj) ~= data.(vname)(nsubj+1:2*nsubj))
        error(['Variable ', vname, ' differs between T0 and T1 rows']);
    end
end

for roi = 1:numel(areas_significant)
    areas_significant{roi}
        if do_lme_correlation == 1

            data_for_corr = table(data.Et_(1:nsubj), ...
                      (data.Sesso(1:nsubj)), ...
                      (data.MCI_s__1_no_0_(1:nsubj)), ...
                      (data.SottotipoMotorio(1:nsubj)), ...
                      data.UMSARSI_1(1:nsubj) - data.UMSARSI(1:nsubj), ...
                      data.DurataMalattia_anni_adjusted(1:nsubj), ...
                      data.(areas_significant{roi})(nsubj+1:end) - data.(areas_significant{roi})(1:nsubj), ...
                      data.(areas_significant{roi})(1:nsubj), ...
                      'VariableNames', {'Et_', 'Sesso', 'MCI_s__1_no_0_', 'SottotipoMotorio', 'deltaUMSARSI','Disease_Duration', ...
                      ['delta_',areas_significant{roi}], ['baseline_',areas_significant{roi}]});

            % mean centering the roi 
            %data.([areas_significant{roi},'_centered']) = zscore(data.(areas_significant{roi})); % - mean(data.(areas_significant{roi}));
            lm_corr_disease = fitlm(data_for_corr, ['Disease_Duration ~ ',['delta_',areas_significant{roi}],' + Et_ + Sesso + MCI_s__1_no_0_ + SottotipoMotorio']);
            lm_corr_umsarsI = fitlm(data_for_corr, ['deltaUMSARSI ~ ',['delta_',areas_significant{roi}],' + Et_ + Sesso + MCI_s__1_no_0_ + SottotipoMotorio']);
            % if any(vifs_lm(lm_corr_disease)>=5)
            %      error(['there is a collinearity issue for model DurataMalattia_anni_adjusted and area ', areas_significant{roi}])
            %  end
            if roi == 1
                recap_table_corr_disease = table({'Intercept','Et_','Sesso','MCI_s__1_no_0_','SottotipoMotorio','ROI'}',  'VariableNames', {'Labels'});
                recap_table_corr_umsarsI = table({'Intercept','Et_','Sesso','MCI_s__1_no_0_','SottotipoMotorio','ROI'}',  'VariableNames', {'Labels'});

            end

        recap_table_corr_disease.([areas_significant{roi},'_estimate']) = round(lm_corr_disease.Coefficients.Estimate,4);
        recap_table_corr_disease.(areas_significant{roi}) = round(lm_corr_disease.Coefficients.pValue,12);

        recap_table_corr_umsarsI.([areas_significant{roi},'_estimate']) = round(lm_corr_umsarsI.Coefficients.Estimate,4);
        recap_table_corr_umsarsI.(areas_significant{roi}) = round(lm_corr_umsarsI.Coefficients.pValue,12);
        end
end

% correlation
close all
fdr_bh(table2array(recap_table_corr_disease(end,3:2:end)),alpha,'pdep','yes')
fdr_bh(table2array(recap_table_corr_umsarsI(end,3:2:end)),alpha,'pdep','yes')

if do_CT == 1 & do_VOL == 0 & fusing_roi == 0
     h_disease = plot_results_corr(recap_table_corr_disease,['CT corr disease duration removed outliers ', num2str(remove_outliers)]);
     set(gca,'XTickLabel',areas_significant)
     h_umsarsI = plot_results_corr(recap_table_corr_umsarsI,['CT corr UMSARSI removed outliers ', num2str(remove_outliers)]);
     set(gca,'XTickLabel',areas_significant)
     exportgraphics(h_disease,['/home/riccardo/codici_progetti/Salerno/distributions_figure/correlations/CT_fusing', num2str(fusing_roi),'_disease_duration_remove_outliers_', num2str(remove_outliers),'.png'], 'Resolution', 300);
     exportgraphics(h_umsarsI,['/home/riccardo/codici_progetti/Salerno/distributions_figure/correlations/CT_fusing', num2str(fusing_roi),'_umsarsI_remove_outliers_', num2str(remove_outliers),'.png'], 'Resolution', 300);

elseif do_CT == 0 & do_VOL == 1
     h_disease = plot_results_corr(recap_table_corr_disease,['VOL corr disease duration removed outliers ', num2str(remove_outliers)]);
     set(gca,'XTickLabel',areas_significant)
     h_umsarsI = plot_results_corr(recap_table_corr_umsarsI,['VOL corr UMSARSI removed outliers ', num2str(remove_outliers)]);
     set(gca,'XTickLabel',areas_significant)
     exportgraphics(h_disease,['/home/riccardo/codici_progetti/Salerno/distributions_figure/correlations/VOL_disease_duration_remove_outliers_', num2str(remove_outliers),'.png'], 'Resolution', 300);
     exportgraphics(h_umsarsI,['/home/riccardo/codici_progetti/Salerno/distributions_figure/correlations/VOL_umsarsI_remove_outliers_', num2str(remove_outliers),'.png'], 'Resolution', 300);

end