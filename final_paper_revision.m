clear
close all
clc

do_CT =0;
do_VOL =1;
atlas_name = 'DKT'; % Schaefer or 'DKT', for cortical thikness stream only 
remove_outliers = 1;
fusing_roi=0;
alpha=0.05;
alpha_bonf=0.0016;
savepath = "/home/riccardo/codici_progetti/MSA_paper/MSA_long_paper/matrices/";

diagnostics = 0;
do_lme_correlation = 1;
NPARC = 400;
NNET = 17;
left_right_merge = true;


if do_VOL == 1 & do_CT == 0
    load('./matrices/tot_stats_vol.mat');
    load('./matrices/tot_stats_brainstem_vol.mat');
    if isequal(tot_stats.Subj_ID,tot_brainstem.Subj_ID) & isequal(tot_stats.SottotipoMotorio,tot_brainstem.SottotipoMotorio)
        data = tot_stats;
        data.Pons = tot_brainstem.Pons;
    end
    areas = {'Putamen', 'CerebellumWhiteMatter','CerebellumCortex','Pons'};
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
            data = [data_schaefer , tot_stats(:,{'Et_','Sesso','SottotipoMotorio','MCI_s__1_no_0_'})];
        else
            error('Subjects IDs are not matching. Please double check')
        end
                
        data.Subj_ID = categorical(data.Subj_ID);
        data.Sesso = categorical(data.Sesso);
        data.Time = categorical(data.Time);
        data.SottotipoMotorio = categorical(data.SottotipoMotorio);
        data.MCI_s__1_no_0_ = categorical(data.MCI_s__1_no_0_);

    else
        error('Atlas name is not valid');
    end
end

lme_list = {'~ Time + MCI_s__1_no_0_ + SottotipoMotorio + Et_ + Sesso + (1|Subj_ID)'}
lme_list_savename={'total_model'};

if remove_outliers == 1
    nsubj_orig = size(data,1);
    [data_clean , outliers_subj] = insert_nan_where_outlier(data, areas); % removing outliers across time points and group (one unique distribution)
    rowsToRemove = contains(string(data.Subj_ID), string(outliers_subj));
    %data(rowsToRemove,:)=[];
    data = data_clean; % for ROI based outlier removal

    %if (nsubj_orig - sum(rowsToRemove)) ~= size(data,1)
    %    error('outliers removal gone wrong')
    %end

    suffix = '_data_clean_outlier_removed';
elseif remove_outliers == 0
    suffix = '';
end

% convert Subj_ID to categorical
data.Subj_ID = categorical(data.Subj_ID);


for type_of_lme = 1:numel(lme_list)
    eta_fenotipo = NaN*ones(numel(areas),1);
    eta_time = NaN*ones(numel(areas),1);
    effect_size_time = NaN*ones(numel(areas),1);
    effect_size_cond =  NaN*ones(numel(areas),1);
    R2_full_model = NaN*ones(numel(areas),1);
    beta_time_standard = NaN*ones(numel(areas),1);
    beta_pheno_standard = NaN*ones(numel(areas),1);

    for roi = 1:numel(areas)
        % loop over type of lme
        areas{roi}
        % remove outlier per ROI 
        data_roi = data;
        badSubj = unique(data_roi.Subj_ID( ismissing(data_roi.(areas{roi}))));
        data_roi( ismember(data_roi.Subj_ID, badSubj), : ) = [];

        lme = fitlme(data_roi, [areas{roi},lme_list{type_of_lme}]);
        % to compute beta relative to baseline 
        % Adjusted baseline mean (population-level)
        mean_bl_ct = table2array(mean(data_roi(data_roi.Time=='T0',areas{roi})));
        mean_bl_marginal_ct = mean(predict(lme, data_roi(data_roi.Time=='T0',:)));
        mean_bl_pheno = table2array(mean(data_roi(data_roi.SottotipoMotorio=='C',areas{roi})));
        mean_bl_marginal_pheno = mean(predict(lme, data_roi(data_roi.SottotipoMotorio=='C',:)));
        if mean_bl_ct - mean_bl_marginal_ct >= 10e-6 | mean_bl_pheno - mean_bl_marginal_pheno >= 10e-6 
            error("the two average at baseline are not equal")
        end
        % redued model without time fixed effect 
        model_string = strrep([areas{roi},lme_list{type_of_lme}], 'Time + ', '  ');  % remove time from the model to calcualte f sqaure
        lme_base = fitlme(data_roi, model_string);
        R2_full =  lme.Rsquared.Ordinary;
        R2_base = lme_base.Rsquared.Ordinary;
        R2_full_model(roi) = R2_full;
        effect_size_time(roi) = (R2_full - R2_base) / (1 - R2_full);
        
        
        if any(vifs(lme)>=5)
            err(['there is a collinearity issue for model UMSARSI and area ', areas{roi}])
        end

        if roi == 1
            recap_table = table( lme.CoefficientNames(1:end)', 'VariableNames', {'Labels'});
        end
        
        time_idx = lme.CoefficientNames == "Time_T1";
        beta_time_standard(roi) = lme.Coefficients.Estimate(time_idx)/mean_bl_ct*100;

        pheno_idx = lme.CoefficientNames == "SottotipoMotorio_P";
        beta_pheno_standard(roi) = lme.Coefficients.Estimate(pheno_idx)/mean_bl_pheno*100;

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

        % effect size for lme 
        res = anova(lme); 

        idx_pheno = find(strcmp(res.Term, 'SottotipoMotorio'), 1);
        F_pheno = res.FStat(idx_pheno);
        df1_pheno = res.DF1(idx_pheno);
        df2_pheno = res.DF2(idx_pheno);
        eta_fenotipo(roi) = (F_pheno*df1_pheno) / ((F_pheno*df1_pheno)+ df2_pheno);

        idx_time   = find(strcmp(res.Term, 'Time'), 1);
        F_time = res.FStat(idx_time);
        df1_time = res.DF1(idx_time);
        df2_time = res.DF2(idx_time);
        eta_time(roi) = (F_time*df1_time) / ((F_time*df1_time)+ df2_time);
    end

    %h_recap=plot_results(recap_table,lme_list{type_of_lme});
    
    disp(['fdr correction for ' , lme_list{type_of_lme}])
    row_time = strcmp(recap_table.Labels, 'Time_T1');
    row_phenotype = strcmp(recap_table.Labels, 'SottotipoMotorio_P');
    
    [h_lme_tot, ~, ~] = fdr_bh(table2array(recap_table(row_time,3:2:end)),alpha,'pdep','yes')
    [h_lme_phenot, ~, ~] = fdr_bh(table2array(recap_table(row_phenotype,3:2:end)),alpha,'pdep','yes')
    if strcmp(atlas_name , 'DKT') == 1
        writetable(recap_table,[ 'recap_MSA_long', suffix, '.tsv'], 'FileType', 'text', 'Delimiter', '\t');
    end
    h_lme_tot_bonf = table2array(recap_table(:,3:2:end))<=alpha_bonf;
    disp('detect if there are more then 1 signficiant values in each roi')
    sum(h_lme_tot_bonf(2:end,:))
    
    % add phenotype and time effect sizes
    recap_table = add_cols(areas, recap_table, 'Phenotype', eta_fenotipo);
    recap_table = add_cols(areas, recap_table, 'Time', eta_time);


if do_CT == 1 && do_VOL == 0
    if strcmp(atlas_name, 'DKT')
        writetable(data, savepath + "CT_MSAlong" + suffix + ".csv");
        writetable(recap_table, savepath + "CT_recap_MSAlong" + suffix + ".tsv", ...
                   'FileType', 'text', 'Delimiter', '\t');
        writetable(recap_table, savepath + "CT_recap_MSAlong" + suffix + ".csv");
        %if remove_outliers == 1
        %    writetable(data_clean, savepath + "CT_MSAlong_data_clean" + suffix + ".csv");
        %end
    elseif strcmp(atlas_name, 'Schaefer')
        writetable(data, savepath + "CT_MSAlong_Schaefer_Parc" + num2str(NPARC) + ...
                         "_NET" + num2str(NNET) + suffix + ".csv");

        writetable(recap_table, savepath + "CT_recap_MSAlong_Schaefer_Parc" + num2str(NPARC) + ...
                              "_NET" + num2str(NNET) + suffix + ".tsv", ...
                   'FileType', 'text', 'Delimiter', '\t');

        writetable(recap_table, savepath + "CT_recap_MSAlong_Schaefer_Parc" + num2str(NPARC) + ...
                              "_NET" + num2str(NNET) + suffix + ".csv");
        %if remove_outliers == 1
        %    writetable(data_clean, savepath + "CT_MSAlong_Schaefer_Parc" + num2str(NPARC) + ...
        %                      "_NET" + num2str(NNET) + "_data_clean" + suffix + ".csv");
        %end
    end

elseif do_CT == 0 && do_VOL == 1
    writetable(data, savepath + "VOL_MSAlong" + suffix + ".csv");
    writetable(recap_table, savepath + "VOL_recap_MSAlong" + suffix + ".tsv", ...
               'FileType', 'text', 'Delimiter', '\t');
    writetable(recap_table, savepath + "VOL_recap_MSAlong" + suffix + ".csv");
    %if remove_outliers == 1
    %        writetable(data_clean, savepath + "VOL_MSAlong_data_clean" + suffix + ".csv");
    %end
end


    
    % vec = 1:numel(areas);
    % figure; hold on; box off;
    % plot(vec, effect_size_marginal,'o')
    % plot(vec(h_lme_tot==1), effect_size_marginal(h_lme_tot==1),'r.')
    % plot(vec , eta_time, 'g*')
    % legend('f square','f sqaure in FDR corrected ROI','partial eta')
    % set(gca,'XTick',vec,'XTickLabel',areas)
end


function recap_table = add_cols(areas, recap_table, name, eta_values)

    for i = 1:numel(areas)
        roi = areas{i};
        eta_colname = [roi '_eta_', name];
        pcol_idx = find(strcmp(recap_table.Properties.VariableNames, roi));
    
        if isempty(pcol_idx)
            error('Column %s not found in recap_table', roi);
        end
    
        eta_col = repmat(eta_values(i), height(recap_table), 1);
    
        % Insert eta column immediately AFTER the p-value column
        recap_table = addvars(recap_table, ...
            eta_col, ...
            'After', pcol_idx, ...
            'NewVariableNames', eta_colname);
    end

end

function R2_marginal = calc_Rsquare(model)

    % Fixed-effects fitted values (X * beta)
    y_fixed = fitted(model,'Conditional',false);
    
    % Residual variance
    sigma2 = model.MSE;
    
    % Random-effect variance (sum of all random-effect variances)
    reVar = sum(model.covarianceParameters{:});
    
    % Variance of fixed-effects prediction
    var_fixed = var(y_fixed,1);
    
    % Marginal R^2
    R2_marginal = var_fixed / (var_fixed + reVar + sigma2);

end

