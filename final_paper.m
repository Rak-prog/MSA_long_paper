clear
close all
clc

do_CT = 0;
do_VOL = 1;
fusing_roi=0;
alpha=0.001;
remove_outliers = 0;
diagnostics = 0;
do_lme_correlation = 1;
% final paper analysis results (it takes the tables results created inlme
% mixed_effect model.m) fitl
% only model total with both interaction MCI and fenotipo
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
lme_list = {'~ Time * MCI_s__1_no_0_ + Time * SottotipoMotorio + Et_ + Sesso + (1|Subj_ID)';...
    '~ Time * UMSARSIcov_centered + SottotipoMotorio + Et_ + Sesso + MCI_s__1_no_0_ + (1|Subj_ID)';...
    '~ Time * disease_duration_centered + SottotipoMotorio + Et_ + Sesso + MCI_s__1_no_0_ + (1|Subj_ID)'}
lme_list_savename={'total_model','umsarsI_model','disease_duration'};

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
            % plotResiduals(lme, 'fitted')
            % plotResiduals(lme, 'probability')
            % %plotDiagnostics(lme, 'leverage')
            disp(lme.Coefficients)
            lme.Coefficients.Estimate ./ lme.Coefficients.SE
        end
    end

    h_recap=plot_results(recap_table,lme_list{type_of_lme});
    if do_CT == 1 & do_VOL == 0
        exportgraphics(h_recap,['/home/riccardo/codici_progetti/Salerno/figure_paper/',lme_list_savename{type_of_lme},...
            'fusing', num2str(fusing_roi), 'remove_outliers_', num2str(remove_outliers),'_lme_ct.png'], 'Resolution', 300);
    elseif do_VOL == 1 & do_CT == 0
        exportgraphics(h_recap,['/home/riccardo/codici_progetti/Salerno/figure_paper/',lme_list_savename{type_of_lme},...
            'fusing', num2str(fusing_roi), 'remove_outliers_', num2str(remove_outliers),'_lme_vol.png'], 'Resolution', 300);
    end
    disp(['fdr correction for ' , lme_list{type_of_lme}])
    [h_lme_tot, ~, ~] = fdr_bh(table2array(recap_table(:,3:2:end)),alpha,'pdep','yes')

    % saving model with fenotipo and MCI interaction i.e. the total modal the paper
    if strcmp('~ Time * MCI_s__1_no_0_ + Time * SottotipoMotorio + Et_ + Sesso + (1|Subj_ID)',lme_list{type_of_lme})...
            && do_CT == 1 && fusing_roi == 1
        writetable(recap_table, ['/home/riccardo/codici_progetti/Salerno/stats_results_lme_final/lme_tot_ct_fused_remove_outliers_', num2str(remove_outliers),'.csv']);
    elseif strcmp('~ Time * MCI_s__1_no_0_ + Time * SottotipoMotorio + Et_ + Sesso + (1|Subj_ID)',lme_list{type_of_lme})...
            && do_CT == 1 && fusing_roi == 0
        writetable(recap_table, ['/home/riccardo/codici_progetti/Salerno/stats_results_lme_final/lme_tot_ct_singleroi_remove_outliers_', num2str(remove_outliers),'.csv']);
    elseif strcmp('~ Time * MCI_s__1_no_0_ + Time * SottotipoMotorio + Et_ + Sesso + (1|Subj_ID)',lme_list{type_of_lme})...
            && do_VOL == 1
        writetable(recap_table, ['/home/riccardo/codici_progetti/Salerno/stats_results_lme_final/lme_tot_vol_remove_outliers_', num2str(remove_outliers),'.csv']);
    end

    % saving UMSARS I
    if strcmp('~ Time * UMSARSIcov_centered + SottotipoMotorio + Et_ + Sesso + MCI_s__1_no_0_ + (1|Subj_ID)',lme_list{type_of_lme})...
            && do_CT == 1 && fusing_roi == 1
        writetable(recap_table, ['/home/riccardo/codici_progetti/Salerno/stats_results_lme_final/lme_tot_ct_fused_umsarsI_remove_outliers_', num2str(remove_outliers),'.csv']);
    elseif strcmp('~ Time * UMSARSIcov_centered + SottotipoMotorio + Et_ + Sesso + MCI_s__1_no_0_ + (1|Subj_ID)',lme_list{type_of_lme})...
            && do_CT == 1 && fusing_roi == 0
        writetable(recap_table, ['/home/riccardo/codici_progetti/Salerno/stats_results_lme_final/lme_tot_ct_singleroi_umsarsI_remove_outliers_', num2str(remove_outliers),'.csv']);
    elseif strcmp('~ Time * UMSARSIcov_centered + SottotipoMotorio + Et_ + Sesso + MCI_s__1_no_0_ + (1|Subj_ID)',lme_list{type_of_lme})...
            && do_VOL == 1
        writetable(recap_table, ['/home/riccardo/codici_progetti/Salerno/stats_results_lme_final/lme_tot_vol_umsarsI_remove_outliers_', num2str(remove_outliers),'.csv']);
    end

    % save disease duration
    if strcmp('~ Time * disease_duration_centered + SottotipoMotorio + Et_ + Sesso + MCI_s__1_no_0_ + (1|Subj_ID)',lme_list{type_of_lme})...
            && do_CT == 1 && fusing_roi == 1
        writetable(recap_table, ['/home/riccardo/codici_progetti/Salerno/stats_results_lme_final/lme_tot_ct_fused_durata_remove_outliers_', num2str(remove_outliers),'.csv']);
    elseif strcmp('~ Time * disease_duration_centered + SottotipoMotorio + Et_ + Sesso + MCI_s__1_no_0_ + (1|Subj_ID)',lme_list{type_of_lme})...
            && do_CT == 1 && fusing_roi == 0
        writetable(recap_table, ['/home/riccardo/codici_progetti/Salerno/stats_results_lme_final/lme_tot_ct_singleroi_durata_remove_outliers_', num2str(remove_outliers),'.csv']);
    elseif strcmp('~ Time * disease_duration_centered + SottotipoMotorio + Et_ + Sesso + MCI_s__1_no_0_ + (1|Subj_ID)',lme_list{type_of_lme})...
            && do_VOL == 1
        writetable(recap_table, ['/home/riccardo/codici_progetti/Salerno/stats_results_lme_final/lme_tot_vol_durata_remove_outliers_', num2str(remove_outliers),'.csv']);
    end
end

% visualize BL vs FUP
data.Subj_ID = string(data.Subj_ID); % converting subj_id column back to cell
dataBL = data(data.Time == 'T0', areas);
dataFUP = data(data.Time == 'T1', areas);
subj_listBL = data(data.Time == 'T0', 'Subj_ID');
subj_listFUP = data(data.Time == 'T1', 'Subj_ID');
min_val = min([ min(min(table2array(dataBL)))   min(min(table2array(dataFUP)))]);
max_val = max([ max(max(table2array(dataBL)))   max(max(table2array(dataFUP)))]);
if do_CT == 0 & do_VOL == 1
    [hvols,~,~,~,~,increasing_subjects]=visualize(table2array(dataBL), table2array(dataFUP), 'BL vs FUP', [min_val max_val], areas,...
        'VOL [%]', table2array(subj_listBL), table2array(subj_listFUP), 1, 'b', 'c')
    exportgraphics(hvols,['/home/riccardo/codici_progetti/Salerno/figure_paper/BLvsFUP_vol_remove_outliers_', num2str(remove_outliers),'.png'], 'Resolution', 300);
elseif do_CT == 1 & do_VOL == 0
    [hct,~,~,~,~,increasing_subjects]=visualize(table2array(dataBL), table2array(dataFUP), 'BL vs FUP', [min_val max_val], areas,...
        'CT [mm]', table2array(subj_listBL), table2array(subj_listFUP), 1, 'b', 'c')
    exportgraphics(hct,['/home/riccardo/codici_progetti/Salerno/figure_paper/BLvsFUP_ct_fusing_roi_',num2str(fusing_roi),'_remove_outliers_',num2str(remove_outliers),'.png'], 'Resolution', 300);
end

% outliers list
max_len = max(cellfun(@numel, increasing_subjects));

% Step 2: Pad shorter cells with empty strings
padded = cell(numel(increasing_subjects), max_len);
for i = 1:numel(increasing_subjects)
    row = increasing_subjects{i};
    padded(i, 1:numel(row)) = row;
end

% Step 3: Convert to table
T = cell2table(padded);
%T.Properties.VariableNames = repmat('', 1, max_len);
T.Properties.RowNames = areas;

% Display the table
if do_CT == 0 & do_VOL == 1
    writetable(T,['/home/riccardo/codici_progetti/Salerno/figure_paper/VOL_subj_increase_remove_outliers_',num2str(remove_outliers),'.csv'],'WriteRowNames', true);
elseif do_CT == 1 & do_VOL == 0
    writetable(T,['/home/riccardo/codici_progetti/Salerno/figure_paper/CT_subj_increase_fusing_',num2str(fusing_roi),'_remove_outliers_',num2str(remove_outliers),'.csv'],'WriteRowNames', true);
end

allCells = table2cell(T);
% Flatten into a single column
flatList = allCells(:);

% Remove empty entries (optional)
flatList = flatList(~cellfun(@isempty, flatList));

% Get unique strings
uniqueStrings = unique(flatList);
clear lme areas
corr_variables = {'disease_duration_centered','UMSARSIcov_centered'};

% add time at FUP to disease duration at FUP
data.DurataMalattia_anni_adjusted(data.Time=='T1') = data.DurataMalattia_anni_(data.Time=='T0') + data.T1_T0_Diff_Days(data.Time=='T0')./365; 
data.DurataMalattia_anni_adjusted(data.Time=='T0') = data.DurataMalattia_anni_(data.Time=='T0'); 
%data.DurataMalattia_anni_adjusted = data.DurataMalattia_anni_adjusted - mean(data.DurataMalattia_anni_adjusted);
%if  ~isequal(round(data.DurataMalattia_anni_adjusted(data.Time=='T1') - data.DurataMalattia_anni_adjusted(data.Time=='T0'),4) , round(data.T1_T0_Diff_Days(data.Time=='T0')./365,4))
%    error('adjusting disease duration is wrong')
%end

% do correaltion for the paper using LME
data.UMSARSIcov_centered = data.UMSARSIcov - mean(data.UMSARSIcov);
for roi = 1:numel(areas_significant)
    areas_significant{roi}
        if do_lme_correlation == 1

            % mean centering the roi 
            data.([areas_significant{roi},'_centered']) = zscore(data.(areas_significant{roi})); % - mean(data.(areas_significant{roi}));
            lme_corr_disease = fitlme(data, ['DurataMalattia_anni_adjusted ~ Time * ', [areas_significant{roi},'_centered'],' + Et_ + Sesso + MCI_s__1_no_0_ + SottotipoMotorio + (1 | Subj_ID)']);
            lme_corr_umsarsI = fitlme(data, ['UMSARSIcov_centered ~ Time * ', [areas_significant{roi},'_centered'],' + Et_ + Sesso + MCI_s__1_no_0_ + SottotipoMotorio + (1 | Subj_ID)']);
            if any(vifs(lme_corr_disease)>=5)
                error(['there is a collinearity issue for model DurataMalattia_anni_adjusted and area ', areas_significant{roi}])
            end
            if roi == 1
                recap_table_corr_disease = table({'Intercept','Et_','Sesso','SottotipoMotorio','MCI_s__1_no_0_','Time_T1','ROI','ROI:Time_T1'}', 'VariableNames', {'Labels'});
                recap_table_corr_umsarsI = table({'Intercept','Et_','Sesso','SottotipoMotorio','MCI_s__1_no_0_','Time_T1','ROI','ROI:Time_T1'}', 'VariableNames', {'Labels'});

            end

        recap_table_corr_disease.([areas_significant{roi},'_estimate']) = round(lme_corr_disease.Coefficients.Estimate,4);
        recap_table_corr_disease.(areas_significant{roi}) = round(lme_corr_disease.Coefficients.pValue,12);

        recap_table_corr_umsarsI.([areas_significant{roi},'_estimate']) = round(lme_corr_umsarsI.Coefficients.Estimate,4);
        recap_table_corr_umsarsI.(areas_significant{roi}) = round(lme_corr_umsarsI.Coefficients.pValue,12);
        end
end

% correlation
close all
fdr_bh(table2array(recap_table_corr_disease(:,3:2:end)),0.001,'pdep','yes')
fdr_bh(table2array(recap_table_corr_umsarsI(:,3:2:end)),0.001,'pdep','yes')

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