% =========================================================================
% Paper: Towards a Digital Twin for the Quantitative Risk Assessment of
%     Domino Effects in Fuel Storage Facilities
% Repository: https://github.com/ramon-quiza/qra-domino-fuel-storage
% File: Escalation model - decides, probabilistically, whether a target
%     tank fails (ignites) due to the thermal radiation it receives
% -------------------------------------------------------------------------
% (c) 2026, Centre for Advanced and Sustainable Manufacturing Studies
%           University of Matanzas
% =========================================================================
% Inputs:
%   I_R_ij - Radiation intensity received by the target tank [kW/m^2]
%   D_ij   - Diameter of the target tank [m]
%   H_ij   - Height of the target tank [m]
% Output:
%   flag   - true if the target tank fails (escalation occurs) in this
%            random trial; false otherwise
% -------------------------------------------------------------------------
% Note on repeated trials: m01 calls this function for every intact tank
% at every escalation step at which new fires have appeared. The trials
% are independent, i.e. a tank that survived an earlier step is tested
% again against the full probability p of the radiation it currently
% receives (see Section 2.4 of the paper for the justification).
% -------------------------------------------------------------------------
% References:
%   [1] Landucci, G., Gubinelli, G., Antonioni, G., Cozzani, V. (2009).
%       The assessment of the damage probability of storage tanks in
%       domino events triggered by fire. Accident Analysis and
%       Prevention, 41(6), 1206-1215.
%       https://doi.org/10.1016/j.aap.2008.05.006
%       (time-to-failure correlation for atmospheric vessels,
%       ln(ttf) = 9.877 - 1.13.ln(I) - 2.667e-5.V, and probit model
%       Y = 9.25 - 1.847.ln(ttf) with ttf in minutes; calibrated so
%       that ttf = 5 min -> ~10% failure probability and
%       ttf = 20 min -> ~90%, accounting for emergency response times)
%   [2] Cozzani, V., Gubinelli, G., Antonioni, G., Spadoni, G.,
%       Zanelli, S. (2005). The assessment of risk caused by domino
%       effect in quantitative area risk analysis. Journal of Hazardous
%       Materials, 127(1-3), 14-30. (probit-based framework for
%       escalation probability in domino effect assessment)
%   [3] Finney, D.J. (1971). Probit Analysis, 3rd ed. Cambridge
%       University Press. (conversion from probit value to probability
%       through the standard normal cumulative distribution function)
% =========================================================================
function flag = escalation(I_R_ij, D_ij, H_ij)
    % Volume of the target (cylindrical) tank [m^3]
    V_ij = pi.*D_ij.^2*H_ij./4;

    % Time to failure [s] of an atmospheric vessel exposed to thermal
    % radiation [1]: decreases with the received radiation intensity
    % and, more weakly, with the vessel volume
    ttf = exp(9.877 - 1.13.*log(I_R_ij) - 2.667e-5.*V_ij);

    % Probit value [1]: the ttf is converted from seconds to minutes
    % (ttf/60); shorter failure times yield higher probit values, i.e.
    % less time available for effective emergency response
    Y = 9.25 - 1.847*log(ttf./60);

    % Failure (escalation) probability [3]: standard normal cumulative
    % distribution function evaluated at (Y - 5), expressed through the
    % error function
    p = 0.5 * (1.0 + erf((Y - 5)/sqrt(2.0)));

    % Bernoulli trial: the tank fails in this scenario with
    % probability p (Monte Carlo sampling of the escalation event)
    flag = (rand(1, 1) < p);