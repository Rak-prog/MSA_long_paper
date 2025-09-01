function [T , subject_outliers]= insert_nan_where_outlier(T, roi_names)
    % Get the table and put NaN where you have an outlier in multiple columns
    % roi_names is a cell array of column names to process
    subject_outliers = {};
    for i = 1:length(roi_names)
        current_col = roi_names{i};
        data = table2array(T(:, current_col));
        med_value = median(data);
        iqr_value = iqr(data)
        
        % Calculate outlier thresholds
        upper_threshold = prctile(data,75) + 1.5*iqr_value;
        lower_threshold = prctile(data,25) - 1.5*iqr_value;
        
        % Find outliers
        idx_outliers = data > upper_threshold | data < lower_threshold;
        
        % Convert table column to appropriate type if needed
        if ~isfloat(T.(current_col))
            T.(current_col) = double(T.(current_col));
        end
        
        % Replace outliers with NaN
        T.(current_col)(idx_outliers) = NaN;
        
        if any(idx_outliers)
            subject_outliers = [subject_outliers ; T.Subj_ID(idx_outliers)];
        end
    end

    subject_outliers = unique(subject_outliers);
end