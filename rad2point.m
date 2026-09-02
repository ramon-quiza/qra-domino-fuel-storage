% =========================================================================
% Paper: Towards a Digital Twin for the Quantitative Risk Assessment of
%     Domino Effects in Fuel Storage Facilities
% Repository: https://github.com/ramon-quiza/qra-domino-fuel-storage
% File: Radiation model - computes the thermal radiation intensity
%     received at an arbitrary point of the terrain from all the
%     burning tanks
% -------------------------------------------------------------------------
% (c) 2026, Centre for Advanced and Sustainable Manufacturing Studies
%           University of Matanzas
% =========================================================================
% Inputs:
%   x_P    - X-axis position of the target point [m]
%   y_P    - Y-axis position of the target point [m]
%   Q_R_i  - External radiation of each tank in the current fire
%            configuration [kW] (zero for intact tanks, > 0 for
%            burning tanks)
%   x      - X-axis centre position of all the tanks [m]
%   y      - Y-axis centre position of all the tanks [m]
%   D      - Diameter of all the tanks [m]
% Output:
%   I_p    - Total radiation intensity received at the point [kW/m^2]
% -------------------------------------------------------------------------
% References:
%   [1] Mudan, K.S. (1984). Thermal radiation hazards from hydrocarbon
%       pool fires. Progress in Energy and Combustion Science, 10(1),
%       59-80. (single point-source radiation model, I = tau_a.Q_R/
%       (4.pi.R^2))
%   [2] CCPS (2000). Guidelines for Chemical Process Quantitative Risk
%       Analysis, 2nd ed. AIChE, New York. (point-source model for
%       consequence analysis; valid for targets farther than ~2.5 pool
%       diameters from the fire)
%   [3] Pietersen, C.M., Huerta, S.M. (1985). Analysis of the LPG
%       incident in San Juan Ixhuatepec, Mexico City. TNO Report,
%       Apeldoorn. (atmospheric transmissivity correlation
%       tau_a = 2.02.(P_w.R)^-0.09; the piecewise fit implemented here
%       matches this form in the intermediate distance range for a
%       fixed water vapour partial pressure)
% =========================================================================
function I_p = rad2point(x_P, y_P, Q_R_i, x, y, D)
I_p = 0;                                         % Accumulated radiation intensity
% -------------------------------------------------------------------------
% Sum the contribution of every burning tank (point-source model)
% -------------------------------------------------------------------------
for t = 1 : length(Q_R_i)
    if Q_R_i(t) > 0                              % Intact tanks contribute nothing
        % Distance from the centre of the emitting tank t to the target
        % point [m], bounded below by the tank radius so that points
        % inside the tank footprint do not yield a singular (R -> 0)
        % or unphysical radiation value
        R = max([norm([x(t) y(t)] - [x_P y_P]), D(t)/2]);
        % Atmospheric transmissivity [3]: fraction of the emitted
        % radiation that reaches the target after absorption by the
        % atmosphere (piecewise empirical correlation on distance,
        % of the form tau_a = a.R^b, cf. Pietersen and Huerta, 1985)
        if (R < 5)
            tau_a = 0.976.*R.^(-0.06);
        elseif (R < 55)
            tau_a = 1.029.*R.^(-0.09);
        else
            tau_a = 1.159.*R.^(-0.12);
        end
        % Point-source model [1,2]: the emitted power spreads over
        % a sphere of radius R, attenuated by the transmissivity
        I_p = I_p + tau_a.*Q_R_i(t)./(4.*pi.*R.^2);
    end
end