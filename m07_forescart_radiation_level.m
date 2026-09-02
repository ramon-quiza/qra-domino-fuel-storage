% =========================================================================
% Paper: Towards a Digital Twin for the Quantitative Risk Assessment of
%     Domino Effects in Fuel Storage Facilities
% Repository: https://github.com/ramon-quiza/qra-domino-fuel-storage
% File: Maps the annual frequency of exceeding a given radiation level
%     at every point of the terrain (exceedance-frequency contour map)
% -------------------------------------------------------------------------
% (c) 2026, Centre for Advanced and Sustainable Manufacturing Studies
%           University of Matanzas
% =========================================================================
% For every point of the computation grid, this script evaluates the
% GRNN surrogate models (m06) to obtain the local Gamma parameters,
% computes the conditional probability that the radiation intensity
% exceeds the level of interest given that an accident has occurred
% (1 - F_Gamma), and multiplies it by the annual accident probability,
% estimated as the fraction n/N of Monte Carlo trials (each equivalent
% to one year of operation, see m01) that resulted in an accident.
% The result is a map of exceedance-frequency contours, in the spirit
% of - but not identical to - the individual risk (fatality) contours
% used in quantitative risk analysis [1]; see Section 4 of the paper.
% -------------------------------------------------------------------------
% References:
%   [1] CCPS (2000). Guidelines for Chemical Process Quantitative Risk
%       Analysis, 2nd ed. AIChE, New York. (individual risk contours as
%       the product of the accident frequency and the conditional
%       probability of the harmful effect)
% =========================================================================
clear all;                                       % Clear workspace variables
clc;                                             % Clear command window

% -------------------------------------------------------------------------
% Tank farm layout: position and geometry of each storage tank
% (needed here only to draw the tank outlines on the map)
% -------------------------------------------------------------------------
x =   [150.0  220.0  177.0  162.5  217.5];       % X-axis centre position [m] (GLOBAL frame)
y =   [150.0  150.0  177.0  240.0  240.0];       % Y-axis centre position [m] (GLOBAL frame)
D =   [ 15.0   10.0   10.0   22.0   22.0];       % Diameter [m]

% -------------------------------------------------------------------------
% Annual accident probability
% -------------------------------------------------------------------------
N = 1e8;                                         % Monte Carlo trials generated in m01 (one per facility-year)
n = 9908;                                        % Trials that resulted in an accident: MUST equal size(id, 1) of the
                                                 % scenarios.mat actually used (n = 10048 in the run reported in the paper)

% rad_level = 2;    % [kW/m^2] - pain threshold for exposed people
% rad_level = 5;    % [kW/m^2] - common safety criterion for people (second-degree burns within tens of seconds of exposure)
% rad_level = 12.5; % [kW/m^2] - piloted ignition of wood and plastics; escalation threshold commonly adopted for domino effect assessment
rad_level = 37.5; % [kW/m^2] - damage to process equipment and steel storage tanks

% Trained GRNN surrogate models (m06): (x, y) -> Gamma parameters
load('./data/grnn_data.mat');
% Grid coordinates X, Y [m] (the Gamma maps A, B are superseded here
% by the network estimates)
load('./data/gamma_data.mat');

% -------------------------------------------------------------------------
% Annual probability of exceeding rad_level at every grid point:
%   P = P(I > rad_level | accident) * P(accident per year)
%     = (1 - F_Gamma(rad_level; a, b)) * n/N                          [1]
% -------------------------------------------------------------------------
for i = 1 : size(X, 1)
    for j = 1 : size(X, 2)
        a = net_a([X(i, j), Y(i, j)].');         % Local Gamma shape parameter [-]
        b = net_b([X(i, j), Y(i, j)].');         % Local Gamma scale parameter [kW/m2]
        p = gamcdf(rad_level, a, b);             % P(I <= rad_level | accident)
        P(i, j) = (1 - p).*n/N;                  % Annual exceedance probability
    end
end

% -------------------------------------------------------------------------
% Plot: exceedance-frequency contours (1e-5 to 9e-5 per year) and tank outlines
% -------------------------------------------------------------------------
fig = figure;
hold on;
contour(X, Y, P, [1 : 9].*1e-5);
% Tank outlines (circles of diameter D centred at (x, y))
theta = linspace(0, 2*pi, 100);
for t = 1 : length(D)
    x_circ = x(t) + D(t)/2.*cos(theta);
    y_circ = y(t) + D(t)/2.*sin(theta);
    plot(x_circ, y_circ, 'k-');
end
axis equal;
colormap('jet');
colorbar;
set(gca, 'XTick', [-50 : 50 : 450]);
set(gca, 'XLim', [-50, 450]);
set(gca, 'YTick', [-50 : 50 : 450]);
set(gca, 'YLim', [-50, 450]);
set(gca, 'CLim', [1, 9].*1e-5);
grid on;
title(sprintf('Radiation = %f kW/m^2', rad_level));

% -------------------------------------------------------------------------
% Save the figure (.fig and .png), labelled by the radiation level
% -------------------------------------------------------------------------
saveas(fig, sprintf('./figures/radlevels/radlevel_%02d.fig', round(rad_level)));
saveas(fig, sprintf('./figures/radlevels/radlevel_%02d.png', round(rad_level)));
