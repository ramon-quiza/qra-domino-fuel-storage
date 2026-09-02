% =========================================================================
% Paper: Towards a Digital Twin for the Quantitative Risk Assessment of
%     Domino Effects in Fuel Storage Facilities
% Repository: https://github.com/ramon-quiza/qra-domino-fuel-storage
% File: Computes the expected number of people affected in a set of
%     buildings, from the probabilistic radiation field of the digital
%     twin, distinguishing sheltered (indoor) and outdoor population
% -------------------------------------------------------------------------
% (c) 2026, Centre for Advanced and Sustainable Manufacturing Studies
%           University of Matanzas
% =========================================================================
% Inputs:
%   x_B    - X-axis position of each building [m]
%   y_B    - Y-axis position of each building [m]
%   P_B    - Usual (time-averaged) number of people in each building
%   net_a  - GRNN surrogate of the Gamma shape parameter (from m06)
%   net_b  - GRNN surrogate of the Gamma scale parameter (from m06)
%   f_out  - Fraction of the residents assumed to be outdoors [-]
%   t_exp  - Effective exposure time of the outdoor population [s]
% Outputs:
%   N_acc  - Expected number of people affected (fatalities) in each
%            building, conditional on the occurrence of an accident
%   p_out  - Probability of harm for an outdoor person at each building
%   p_in   - Probability of harm for an indoor person at each building
% -------------------------------------------------------------------------
% Vulnerability model:
%   - OUTDOOR population: probit of the thermal radiation dose,
%     V = t_exp.I^(4/3) (I in W/m2), Y = -14.9 + 2.56.ln(V/1e4),
%     following Eisenberg et al. [1] as adopted by the TNO Green
%     Book [2]. The probability of harm is obtained by averaging the
%     probit response over the local Gamma distribution of the
%     radiation intensity provided by the digital twin.
%   - INDOOR (sheltered) population: buildings shield their occupants
%     from direct radiation, so harm is driven by the ignition of the
%     building itself; following the convention of the Dutch QRA
%     practice [3], full harm is assumed when the outdoor radiation
%     exceeds I_ign = 35 kW/m2 and no harm below it, so the probability
%     of harm reduces to the Gamma exceedance of that threshold.
%   - The effective exposure time of people outdoors accounts for the
%     time needed to reach shelter or escape; values of 20-30 s are
%     recommended in the literature [2, 3], with t_exp = 20 s as the
%     usual default.
% -------------------------------------------------------------------------
% References:
%   [1] Eisenberg, N.A., Lynch, C.J., Breeding, R.J. (1975).
%       Vulnerability Model: A Simulation System for Assessing Damage
%       Resulting from Marine Spills. Report CG-D-136-75, US Coast
%       Guard, Washington DC. (thermal dose probit for people)
%   [2] TNO (1992). Methods for the Determination of Possible Damage
%       to People and Objects Resulting from Releases of Hazardous
%       Materials (Green Book), CPR 16E. Committee for the Prevention
%       of Disasters, The Hague. (dose-based vulnerability of people;
%       effective exposure duration)
%   [3] RIVM (2009). Reference Manual Bevi Risk Assessments, version
%       3.2. National Institute of Public Health and the Environment,
%       Bilthoven. (protection offered by buildings; 35 kW/m2 as the
%       threshold above which indoor lethality is assumed)
% =========================================================================
function [N_acc, p_out, p_in] = people_affected(x_B, y_B, P_B, ...
                                                net_a, net_b, f_out, t_exp)

    I_ign = 35;                                  % Radiation threshold for building ignition [kW/m2] [3]

    % Preallocate the per-building outputs
    N_acc = zeros(size(P_B));
    p_out = zeros(size(P_B));
    p_in  = zeros(size(P_B));

    for k = 1 : length(P_B)
        % Local Gamma parameters of the radiation intensity at the
        % building, provided by the GRNN surrogates of the digital twin
        a = net_a([x_B(k), y_B(k)].');           % Shape parameter [-]
        b = net_b([x_B(k), y_B(k)].');           % Scale parameter [kW/m2]

        % -----------------------------------------------------------------
        % Outdoor population: probability of harm averaged over the
        % local distribution of the radiation intensity,
        %   p_out = Int_0^Inf Phi(Y(I) - 5) f_Gamma(I; a, b) dI
        % with the thermal dose probit Y(I) [1, 2] (I converted from
        % kW/m2 to W/m2 inside the dose)
        % -----------------------------------------------------------------
        probit  = @(I) -14.9 + 2.56.*log(t_exp.*(1e3.*I).^(4/3)./1e4);
        harm    = @(I) normcdf(probit(I) - 5);   % Probit -> probability [1]
        p_out(k) = integral(@(I) harm(I).*gampdf(I, a, b), 0, Inf);

        % -----------------------------------------------------------------
        % Indoor (sheltered) population: harm only if the radiation
        % exceeds the building ignition threshold [3],
        %   p_in = 1 - F_Gamma(I_ign; a, b)
        % -----------------------------------------------------------------
        p_in(k) = 1 - gamcdf(I_ign, a, b);

        % -----------------------------------------------------------------
        % Expected number of people affected in the building,
        % conditional on the occurrence of an accident
        % -----------------------------------------------------------------
        N_acc(k) = P_B(k).*(f_out.*p_out(k) + (1 - f_out).*p_in(k));
    end
end