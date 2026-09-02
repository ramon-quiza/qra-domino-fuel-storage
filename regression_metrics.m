% =========================================================================
% Paper: Towards a Digital Twin for the Quantitative Risk Assessment of
%     Domino Effects in Fuel Storage Facilities
% Repository: https://github.com/ramon-quiza/qra-domino-fuel-storage
% File: Computes the goodness-of-fit metrics (R2, RMSE and maximum
%     absolute error) of a surrogate model, e.g. the GRNN networks of
%     the digital twin, from the observed and predicted values
% -------------------------------------------------------------------------
% (c) 2026, Centre for Advanced and Sustainable Manufacturing Studies
%           University of Matanzas
% =========================================================================
% Inputs:
%   y      - Observed (reference) values, e.g. the Gamma parameters
%            fitted by maximum likelihood in m05
%   y_hat  - Predicted values, e.g. the GRNN estimates at the same
%            points (same size and ordering as y)
% Outputs:
%   R2     - Coefficient of determination [-]:
%            R2 = 1 - SS_res/SS_tot, where SS_res is the sum of squared
%            residuals and SS_tot the total sum of squares about the
%            mean of the observations; R2 = 1 indicates a perfect
%            reproduction, R2 = 0 a model no better than the mean, and
%            R2 < 0 a model worse than the mean [1]
%   RMSE   - Root mean square error, in the units of y
%   e_max  - Maximum absolute error, in the units of y
% -------------------------------------------------------------------------
% Notes:
%   - Both inputs are flattened internally, so vectors, grids (e.g. the
%     101 x 101 parameter maps) or any equally-sized arrays can be
%     passed directly.
%   - Pairs containing NaN in either input are excluded from the
%     computation (e.g. validation points masked out of the training).
% -------------------------------------------------------------------------
% References:
%   [1] Montgomery, D.C., Runger, G.C. (2014). Applied Statistics and
%       Probability for Engineers, 6th ed. Wiley, Hoboken. (definition
%       of the coefficient of determination and regression error
%       metrics)
% =========================================================================
function [R2, RMSE, e_max] = regression_metrics(y, y_hat)

    % ---------------------------------------------------------------------
    % Input conditioning: flatten to column vectors and check sizes
    % ---------------------------------------------------------------------
    y     = y(:);
    y_hat = y_hat(:);
    if length(y) ~= length(y_hat)
        error('regression_metrics:size', ...
              'Observed and predicted arrays must have the same number of elements.');
    end

    % Exclude pairs with NaN in either input (e.g. masked points)
    valid = ~isnan(y) & ~isnan(y_hat);
    y     = y(valid);
    y_hat = y_hat(valid);

    % ---------------------------------------------------------------------
    % Metrics
    % ---------------------------------------------------------------------
    res    = y - y_hat;                          % Residuals
    SS_res = sum(res.^2);                        % Sum of squared residuals
    SS_tot = sum((y - mean(y)).^2);              % Total sum of squares

    R2    = 1 - SS_res./SS_tot;                  % Coefficient of determination [-] [1]
    RMSE  = sqrt(mean(res.^2));                  % Root mean square error
    e_max = max(abs(res));                       % Maximum absolute error
end
