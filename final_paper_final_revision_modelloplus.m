clear
close all
clc

do_CT = 0;
do_VOL = 1;
fusing_roi=0;
alpha=0.05;
alpha_bonf=0.0016;
remove_outliers = 0;
diagnostics = 1;
do_lme_correlation = 1;
% final paper analysis results (it takes the tables results created inlme
% mixed_effect model.m) fitl
% only model total with both interacti.on MCI and fenotipo
% UMSARSI only
% and disease duration (interesting for Pellecchia)

if do_VOL == 1 & do_CT == 0
    load('./matrices/tot_stats_vol.mat');
    load('./matrices/tot_stats_brainstem_vol.mat');
    if isequal(tot_stats.Subj_ID,tot_brainstem.Subj_ID) & isequal(tot_stats.SottotipoMotorio,tot_brainstem.SottotipoMotorio)
        data = tot_stats;
        data.Pons = tot_brainstem.Pons;
    end
    areas = {'Putamen', 'CerebellumWhiteMatter','CerebellumCortex','Pons'};
    areas_significant = {'Putamen','Pons'}; % the areas that survived correction for multiple comparison for the volume at p < 0.001
elseif   do_VOL ==  0 & do_CT == 1
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
        areas_significant = {'postcentral','posteriorcingulate','isthmuscingulate'};
    end
end
% Time * MCI_s__1_no_0_ + Time * SottotipoMotorio + Et_ + Sesso + (1|Subj_ID)
lme_list = {'~ Time + MCI_s__1_no_0_ + SottotipoMotorio + Et_ + Sesso + (1|Subj_ID)'}
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
            plotResiduals(lme, 'fitted')

            figure;
            qqplot(res_std);
            title(['LME: Q-Q plot of standardized residuals', areas{roi}]);
            %plotResiduals(lme, 'probability')

            res_std_centered = (res_std - mean(res_std)) / std(res_std);
            [h_ks, p_ks] = kstest(res_std_centered);
            disp([h_ks, p_ks]);
            if h_ks ~= 0
                disp(['KS TEST IS SIGNFICANT FOR ROI ', areas{roi}])
            end

            F = fitted(lme);
            R = response(lme);
            figure();
            plot(R,F,'rx')
            xlabel('Response')
            ylabel('Fitted')
            title(['LME: Observed Response vs Fitted Response plot ', areas{roi}]);

            close all;
        end
    end

    %h_recap=plot_results(recap_table,lme_list{type_of_lme});
    
    disp(['fdr correction for ' , lme_list{type_of_lme}])
    row_group = strcmp(recap_table.Labels, 'Time_T1');
    [h_lme_tot, ~, ~] = fdr_bh(table2array(recap_table(row_group,3:2:end)),alpha,'pdep','yes')
    writetable(recap_table,[ 'recap_MSA_long.tsv'], 'FileType', 'text', 'Delimiter', '\t');

    h_lme_tot_bonf = table2array(recap_table(:,3:2:end))<=alpha_bonf;
    disp('detect if there are more then 1 signficiant values in each roi')
    sum(h_lme_tot_bonf(2:end,:))
    % % saving model with fenotipo and MCI interaction i.e. the total modal the paper
    % if strcmp('~ Time * MCI_s__1_no_0_ + Time * SottotipoMotorio + Et_ + Sesso + (1|Subj_ID)',lme_list{type_of_lme})...
    %         && do_CT == 1 && fusing_roi == 1
    %     writetable(recap_table, ['/home/riccardo/codici_progetti/Salerno/stats_results_lme_final/lme_tot_ct_fused_remove_outliers_', num2str(remove_outliers),'.csv']);
    % elseif strcmp('~ Time * MCI_s__1_no_0_ + Time * SottotipoMotorio + Et_ + Sesso + (1|Subj_ID)',lme_list{type_of_lme})...
    %         && do_CT == 1 && fusing_roi == 0
    %     writetable(recap_table, ['/home/riccardo/codici_progetti/Salerno/stats_results_lme_final/lme_tot_ct_singleroi_remove_outliers_', num2str(remove_outliers),'.csv']);
    % elseif strcmp('~ Time * MCI_s__1_no_0_ + Time * SottotipoMotorio + Et_ + Sesso + (1|Subj_ID)',lme_list{type_of_lme})...
    %         && do_VOL == 1
    %     writetable(recap_table, ['/home/riccardo/codici_progetti/Salerno/stats_results_lme_final/lme_tot_vol_remove_outliers_', num2str(remove_outliers),'.csv']);
    % end

    % % saving UMSARS I
    % if strcmp('~ Time * UMSARSIcov_centered + SottotipoMotorio + Et_ + Sesso + MCI_s__1_no_0_ + (1|Subj_ID)',lme_list{type_of_lme})...
    %         && do_CT == 1 && fusing_roi == 1
    %     writetable(recap_table, ['/home/riccardo/codici_progetti/Salerno/stats_results_lme_final/lme_tot_ct_fused_umsarsI_remove_outliers_', num2str(remove_outliers),'.csv']);
    % elseif strcmp('~ Time * UMSARSIcov_centered + SottotipoMotorio + Et_ + Sesso + MCI_s__1_no_0_ + (1|Subj_ID)',lme_list{type_of_lme})...
    %         && do_CT == 1 && fusing_roi == 0
    %     writetable(recap_table, ['/home/riccardo/codici_progetti/Salerno/stats_results_lme_final/lme_tot_ct_singleroi_umsarsI_remove_outliers_', num2str(remove_outliers),'.csv']);
    % elseif strcmp('~ Time * UMSARSIcov_centered + SottotipoMotorio + Et_ + Sesso + MCI_s__1_no_0_ + (1|Subj_ID)',lme_list{type_of_lme})...
    %         && do_VOL == 1
    %     writetable(recap_table, ['/home/riccardo/codici_progetti/Salerno/stats_results_lme_final/lme_tot_vol_umsarsI_remove_outliers_', num2str(remove_outliers),'.csv']);
    % end

    % % save disease duration
    % if strcmp('~ Time * disease_duration_centered + SottotipoMotorio + Et_ + Sesso + MCI_s__1_no_0_ + (1|Subj_ID)',lme_list{type_of_lme})...
    %         && do_CT == 1 && fusing_roi == 1
    %     writetable(recap_table, ['/home/riccardo/codici_progetti/Salerno/stats_results_lme_final/lme_tot_ct_fused_durata_remove_outliers_', num2str(remove_outliers),'.csv']);
    % elseif strcmp('~ Time * disease_duration_centered + SottotipoMotorio + Et_ + Sesso + MCI_s__1_no_0_ + (1|Subj_ID)',lme_list{type_of_lme})...
    %         && do_CT == 1 && fusing_roi == 0
    %     writetable(recap_table, ['/home/riccardo/codici_progetti/Salerno/stats_results_lme_final/lme_tot_ct_singleroi_durata_remove_outliers_', num2str(remove_outliers),'.csv']);
    % elseif strcmp('~ Time * disease_duration_centered + SottotipoMotorio + Et_ + Sesso + MCI_s__1_no_0_ + (1|Subj_ID)',lme_list{type_of_lme})...
    %         && do_VOL == 1
    %     writetable(recap_table, ['/home/riccardo/codici_progetti/Salerno/stats_results_lme_final/lme_tot_vol_durata_remove_outliers_', num2str(remove_outliers),'.csv']);
    % end
end
