function vif_tbl = vifs_lm(lm)
%VIFS_LM Compute variance inflation factors for a fitlm LinearModel.
%
%   vif_tbl = vifs_lm(lm)
%
%   lm      : LinearModel object returned by fitlm
%   vif_tbl : table with predictor names and VIF values
%
%   Notes:
%   - Intercept is excluded.
%   - Works with continuous and dummy-coded predictors.

    % Design matrix and coefficient names
    X = lm.DesignMatrix;               % N x p
    names = lm.CoefficientNames(:);    % p x 1 cell

    % Drop intercept column
    isIntercept = strcmp(names,'(Intercept)');
    X = X(:, ~isIntercept);
    names = names(~isIntercept);

    nPred = numel(names);
    vifs = nan(nPred,1);

    for j = 1:nPred
        y = X(:, j);                             % predictor to regress on others
        Xj = X(:, setdiff(1:nPred, j));          % all other predictors

        % add intercept for the auxiliary regression
        Xj = [ones(size(Xj,1),1), Xj];

        % ordinary least squares
        beta = Xj \ y;
        yhat = Xj * beta;

        % R^2 of this auxiliary regression
        RSS = sum((y - yhat).^2);
        TSS = sum((y - mean(y)).^2);
        R2  = 1 - RSS / TSS;

        vifs(j) = 1 / (1 - R2);
    end

    vif_tbl = table(string(names), vifs, ...
                    'VariableNames', {'Predictor','VIF'});
end
