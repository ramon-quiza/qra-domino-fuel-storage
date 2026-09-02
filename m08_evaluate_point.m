% =========================================================================
% Paper: Towards a Digital Twin for the Quantitative Risk Assessment of
%     Domino Effects in Fuel Storage Facilities
% Repository: https://github.com/ramon-quiza/qra-domino-fuel-storage
% File: Evaluates the annual exceedance probability of several
%     radiation levels at an arbitrary point of the terrain
% -------------------------------------------------------------------------
% (c) 2026, Centre for Advanced and Sustainable Manufacturing Studies
%           University of Matanzas
% =========================================================================
% This script illustrates the use of the digital twin at the level of a
% single receptor: given an arbitrary point (x, y) of the terrain (e.g.
% the location of a building or a workplace), the GRNN surrogate models
% (m06) provide the local Gamma parameters of the radiation intensity,
% from which the annual probability of exceeding a set of radiation
% thresholds is computed as in m07 (conditional exceedance probability
% times the annual accident probability n/N) [1].
% -------------------------------------------------------------------------
% References:
%   [1] CCPS (2000). Guidelines for Chemical Process Quantitative Risk
%       Analysis, 2nd ed. AIChE, New York. (individual risk at a given
%       location as the product of the accident frequency and the
%       conditional probability of the harmful effect)
% =========================================================================
clear all;                                       % Clear workspace variables
clc;                                             % Clear command window

% -------------------------------------------------------------------------
% Tank farm layout (kept for reference; only the surrogate models are
% actually needed to evaluate a point)
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

rad_level = 2;  % [kW/m^2] (not used; the thresholds are given in R below)

% Trained GRNN surrogate models (m06): (x, y) -> Gamma parameters
load('./data/grnn_data.mat');

% -------------------------------------------------------------------------
% Receptor of interest [m], in the GLOBAL frame (here, the position of
% the administrative office building B1 of the complex; note that this
% reuses -shadows- the variable names x, y, no longer needed as tank
% positions)
% -------------------------------------------------------------------------
x = 215.3;
y = 304.0;

% Local Gamma parameters of the radiation intensity at the receptor
a = net_a([x, y].');                             % Shape parameter [-]
b = net_b([x, y].');                             % Scale parameter [kW/m2]

% -------------------------------------------------------------------------
% Annual probability of exceeding each radiation threshold at the
% receptor: f = (1 - F_Gamma(R; a, b)) * n/N                          [1]
% -------------------------------------------------------------------------
R = [2 5 12.5 37.5];                             % Radiation thresholds of interest [kW/m2]
f = (1 - gamcdf(R, a, b)).*n/N;                  % Annual exceedance probability of each threshold

% Plot: exceedance probability against the radiation threshold
plot(R, f);
