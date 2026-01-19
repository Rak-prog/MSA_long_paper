function vifs = vifs(mdl)
    % Get the design matrix (may contain NaNs)
    X = mdl.designMatrix('Fixed');
    
    % Check if first column is an intercept (all 1s)
    if all(X(:,1) == 1)
        X = X(:, 2:end);  % Remove intercept
    end
    
    % Initialize VIFs with NaN (in case computation fails for some columns)
    vifs = NaN(size(X, 2), 1);
    
    % Compute VIF for each predictor (handling NaNs column-wise)
    for i = 1:size(X, 2)
        % Skip if column is all NaN or constant (avoid division by zero)
        if all(isnan(X(:,i))) || var(X(:,i), 'omitnan') == 0
            continue;
        end
        
        % Find complete cases (rows without NaN in current predictor & others)
        valid_rows = ~isnan(X(:,i));
        X_clean = X(valid_rows, :);
        
        % Check if enough samples remain (minimum 3 for correlation)
        if size(X_clean, 1) < 3
            warning('Not enough non-NaN samples for column %d. Skipping.', i);
            continue;
        end
        
        % Compute correlation matrix (omit NaN pairs)
        R = corr(X_clean, 'rows', 'complete');
        
        % Check if R is invertible
        if rcond(R) < eps(class(R))
            warning('Correlation matrix near-singular for column %d. Skipping.', i);
            continue;
        end
        
        % Calculate VIF for this predictor
        Rinv = inv(R);
        vifs(i) = Rinv(i,i);
    end
end