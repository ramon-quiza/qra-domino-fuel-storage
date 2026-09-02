% =========================================================================
% Paper: Towards a Digital Twin for the Quantitative Risk Assessment of
%     Domino Effects in Fuel Storage Facilities
% Repository: https://github.com/ramon-quiza/qra-domino-fuel-storage
% File: Pool fire model - computes the heat radiated by a burning tank
% -------------------------------------------------------------------------
% (c) 2026, Centre for Advanced and Sustainable Manufacturing Studies
%           University of Matanzas
% =========================================================================
% Inputs:
%   D    - Tank (pool) diameter [m]
%   T_b  - Fuel boiling point [oC]
%   H_c  - Fuel heat of combustion [kJ/kg]
%   H_v  - Fuel heat of vaporization [kJ/kg]
%   c_P  - Fuel specific heat [kJ/(kg.K)]
% Outputs:
%   Q_R  - External (radiated) heat release rate of the pool fire [kW]
%   m    - Mass burning rate per unit area [kg/(m^2.s)]
%   H_f  - Flame height [m] (informative only: it is not used by the
%          point-source radiation model of the paper)
% -------------------------------------------------------------------------
% Model parameters (Section 2.2 of the paper):
%   k        = 8.76e-4 kg/(m^2.s)  burning-rate constant [1]
%   T_0      = 20 oC               ambient temperature
%   lambda_c = 0.92                combustion efficiency
%   chi_r    = 0.30                flame radiative fraction, held
%                                  constant (a conservative choice for
%                                  large pools; see the paper)
%   delta    = 5 m                 enlargement of the burning pool
%                                  beyond the tank shell (D + 5)
% -------------------------------------------------------------------------
% References:
%   [1] Chen et al. (2023). Pool fire burning characteristics and risks
%       under wind-free conditions: State-of-the-art. Fire Safety
%       Journal, 136, 103755. (reports the fitting coefficient
%       k = 8.76e-4 kg/(m^2.s))
%       https://doi.org/10.1016/j.firesaf.2023.103755
%   [2] Zabetakis, M.G., Burgess, D.S. (1961). Research on the Hazards
%       Associated with the Production and Handling of Liquid Hydrogen.
%       Report BM-RI-5707, US Bureau of Mines, Washington DC. (original
%       energy-balance form of the burning-rate correlation)
%   [3] Thomas, P.H. (1963). The size of flames from natural fires.
%       Symposium (International) on Combustion, 9(1), 844-859.
%       (flame-height correlation H_f/D = 42.[m''/(rho_a.sqrt(g.D))]^0.61)
%   [4] Mudan, K.S. (1984). Thermal radiation hazards from hydrocarbon
%       pool fires. Progress in Energy and Combustion Science, 10(1),
%       59-80. (radiative fraction chi_r and point-source heat release)
% =========================================================================
function [Q_R, m, H_f] = poolfire(D, T_b, H_c, H_v, c_P)
    % ---------------------------------------------------------------------
    % Physical constants and model parameters
    % ---------------------------------------------------------------------
    g = 9.81;             % Acceleration due to gravity [m/s^2]
    C_d = 0.62;           % Discharge coefficient (not used by this model; kept for reference)
    rho_a = 1.225;        % Air mass density [kg/m^3]
    k = 8.76e-4;          % Constant [kg/(m^2.s)]
    T_0 = 20;             % Air temperature [oC]
    lambda_c = 0.92;      % Combustion efficiency
    chi_r = 0.30;         % Flame radiative fraction

    % ---------------------------------------------------------------------
    % Pool fire model
    % ---------------------------------------------------------------------
    % Mass burning rate per unit area [kg/(m^2.s)] [1,2]: fuel evaporation rate
    % driven by the ratio between the combustion heat and the energy
    % needed to heat the fuel from ambient temperature to its boiling
    % point and vaporize it
    m = k.*H_c./(H_v + c_P.*(T_b - T_0));

    % Flame height [m] - Thomas correlation [3], based on the
    % dimensionless burning rate (Froude-number scaling). Informative
    % only: the point-source model does not use the flame geometry
    H_f = 42.*D.*(m./(rho_a.*sqrt(g.*D))).^0.61;

    % Radiated power [kW] [4]: fraction chi_r of the total combustion
    % power emitted as thermal radiation, over the pool area enlarged by
    % a delta = 5 m margin around the tank shell (Eq. (3) of the paper)
    Q_R = chi_r.*lambda_c.*m.*H_c.*pi.*(D + 5).^2./4;