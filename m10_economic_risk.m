% =========================================================================
% Paper: Towards a Digital Twin for the Quantitative Risk Assessment of
%     Domino Effects in Fuel Storage Facilities
% Repository: https://github.com/ramon-quiza/qra-domino-fuel-storage
% File: Computes the Expected Annual Loss (EAL) of every asset of the
%     complex from the probabilistic radiation field of the digital twin
% -------------------------------------------------------------------------
% (c) 2026, Centre for Advanced and Sustainable Manufacturing Studies
%           University of Matanzas
% =========================================================================
% For every asset (storage tanks and their fuel inventory, industrial
% machinery and buildings), this script evaluates the GRNN surrogate
% models (m06) to obtain the local Gamma distribution of the radiation
% intensity, computes the expected damage fraction through a stepped
% vulnerability function D(I) defined per asset class, and obtains the
% Expected Annual Loss as
%
%   EAL_i = (n/N) . V_i . E[D_i],
%   E[D_i] = Sum_k d_k . [F_Gamma(I_k+1) - F_Gamma(I_k)],
%
% i.e. the annual accident probability times the asset value times the
% expected damage fraction, the latter reduced to differences of the
% Gamma cumulative distribution function evaluated at the damage
% thresholds [1]. All positions are given in the GLOBAL coordinate
% system of the complex (see the layout figure and data tables), the
% same frame in which the whole pipeline (m01-m09) and the surrogate
% models operate (OFFSET = 0).
% -------------------------------------------------------------------------
% Damage thresholds of the stepped vulnerability functions [kW/m2]:
%   - Robust equipment and storage tanks: partial damage (30%) between
%     12.5 (piloted ignition of ancillary materials, seals and
%     coatings) and 37.5 (damage to process equipment and steel tanks);
%     total loss above 37.5 [2, 3].
%   - Sensitive equipment (electrical, instrumented): partial damage
%     (30%) between 5 and 12.5; total loss above 12.5, owing to the
%     vulnerability of electronics, cabling and elastomers [2, 3].
%   - Buildings: light damage (10%, glazing and cladding) between 5
%     and 12.5; total loss above 12.5, a conservative assumption
%     justified by the hours-long duration of tank fires, which makes
%     piloted ignition of the envelope the governing mechanism [2, 3].
%   - Fuel inventory: lost with the tank (total loss above 37.5).
% These fractions are engineering assumptions of this work and can be
% refined; the thresholds are the standard consequence-analysis values.
% -------------------------------------------------------------------------
% References:
%   [1] CCPS (2000). Guidelines for Chemical Process Quantitative Risk
%       Analysis, 2nd ed. AIChE, New York. (risk as frequency times
%       expected consequence; damage thresholds for equipment)
%   [2] Raj, P.K. (2008). A review of the criteria for people exposure
%       to radiant heat flux from fires. Journal of Hazardous
%       Materials, 159(1), 61-71. (harm criteria for thermal radiation)
%   [3] Mudan, K.S. (1984). Thermal radiation hazards from hydrocarbon
%       pool fires. Progress in Energy and Combustion Science, 10(1),
%       59-80. (thermal radiation damage to structures and equipment)
% =========================================================================
clear all;                                       % Clear workspace variables
clc;                                             % Clear command window

% -------------------------------------------------------------------------
% Coordinate systems
% -------------------------------------------------------------------------
OFFSET = 0;                                      % Frame offset [m] between the asset coordinates and the surrogate frame (0: same global frame)

% -------------------------------------------------------------------------
% Annual accident probability (see m01 and m07)
% -------------------------------------------------------------------------
N = 1e8;                                         % Monte Carlo trials generated in m01 (one per facility-year)
n = 9908;                                        % Trials that resulted in an accident: MUST equal size(id, 1) of the
                                                 % scenarios.mat actually used (n = 10048 in the run reported in the paper)

% -------------------------------------------------------------------------
% Tank farm layout (GLOBAL coordinates) and fuel inventory model
% -------------------------------------------------------------------------
x_T = [150.0  220.0  177.0  162.5  217.5];       % X-axis centre position [m]
y_T = [150.0  150.0  177.0  240.0  240.0];       % Y-axis centre position [m]
D_T = [ 15.0   10.0   10.0   22.0   22.0];       % Diameter [m]
H_T = [ 12.0    8.0   10.0   12.0   12.0];       % Height [m]
V_T = [0.70   0.35   0.35   1.19   1.19].*1e6;   % Replacement value of the tank [EUR]
c_F = [ 500    700    700    700    700];        % Fuel unit value [EUR/m3] (naphtha, jet, jet, diesel, diesel)
FILL = 0.80;                                     % Typical fill fraction of the tanks [-]
% Fuel inventory value of each tank [EUR]
V_F = c_F.*FILL.*pi.*D_T.^2.*H_T./4;

% -------------------------------------------------------------------------
% Asset inventory (GLOBAL coordinates)
% Classes of the stepped vulnerability functions:
%   1 - robust equipment and storage tanks
%   2 - sensitive (electrical/instrumented) equipment
%   3 - buildings
%   4 - fuel inventory (lost with the tank)
% -------------------------------------------------------------------------
assets = {
%    Id      x [m]   y [m]   Value [EUR]  Class
    'T1',    150.0,  150.0,  V_T(1),      1;
    'T2',    220.0,  150.0,  V_T(2),      1;
    'T3',    177.0,  177.0,  V_T(3),      1;
    'T4',    162.5,  240.0,  V_T(4),      1;
    'T5',    217.5,  240.0,  V_T(5),      1;
    'F1',    150.0,  150.0,  V_F(1),      4;     % Fuel inventory of T1
    'F2',    220.0,  150.0,  V_F(2),      4;     % Fuel inventory of T2
    'F3',    177.0,  177.0,  V_F(3),      4;     % Fuel inventory of T3
    'F4',    162.5,  240.0,  V_F(4),      4;     % Fuel inventory of T4
    'F5',    217.5,  240.0,  V_F(5),      4;     % Fuel inventory of T5
    'M1',    145.6,  195.9,  0.35e6,      1;     % Product pumping and transfer station
    'M2',    237.9,  195.9,  0.35e6,      1;     % Product pumping and transfer station
    'M3',    145.6,  268.6,  1.50e6,      2;     % Road tanker loading gantry (metering/automation)
    'M4',    237.9,  268.6,  0.55e6,      2;     % Electrical substation and motor control centre
    'M5',    209.1,  374.2,  0.30e6,      1;     % Goliath crane
    'M6',    239.0,  374.2,  0.15e6,      2;     % Compressed air plant
    'B1',    215.3,  304.0,  5.40e6,      3;     % Administrative office building
    'B2',    159.1,  374.2, 32.90e6,      3;     % Heavy engineering works
    'B3',    292.0,  369.5,  0.56e6,      3;     % 2-storey detached house
    'B4',    317.0,  369.5,  0.56e6,      3;     % 2-storey detached house
    'B5',    342.0,  369.5,  0.56e6,      3;     % 2-storey detached house
    'B6',    367.0,  369.5,  0.56e6,      3;     % 2-storey detached house
    'B7',     96.9,  369.5,  0.56e6,      3;     % 2-storey detached house
    'B8',     71.9,  369.5,  0.56e6,      3;     % 2-storey detached house
    'B9',    327.5,  276.6, 18.00e6,      3;     % Warehouse
    'B10',   337.5,  204.6, 22.75e6,      3;     % Warehouse
    'B11',   337.5,  159.6, 22.75e6,      3;     % Warehouse
    'B12',    83.7,  250.0, 22.75e6,      3;     % Warehouse
    'B13',   338.8,   80.6, 13.20e6,      3;     % 5-storey block of flats
    'B14',   338.8,   45.6, 13.20e6,      3;     % 5-storey block of flats
    'B15',   338.8,   10.6, 13.20e6,      3;     % 5-storey block of flats
    'B16',   248.8,   80.6, 13.20e6,      3;     % 5-storey block of flats
    'B17',   248.8,   45.6, 13.20e6,      3;     % 5-storey block of flats
    'B18',   248.8,   10.6, 13.20e6,      3;     % 5-storey block of flats
    'B19',   158.8,   80.6, 13.20e6,      3;     % 5-storey block of flats
    'B20',   158.8,   45.6, 13.20e6,      3;     % 5-storey block of flats
    'B21',   158.8,   10.6, 13.20e6,      3;     % 5-storey block of flats
    'B22',    66.0,   21.5, 13.20e6,      3;     % 5-storey block of flats
};

% -------------------------------------------------------------------------
% Stepped vulnerability functions: bands{c} = [I_low, I_high, fraction]
% (thresholds in kW/m2; fractions of the asset value lost)
% -------------------------------------------------------------------------
bands = cell(4, 1);
bands{1} = [12.5   37.5  0.30;                   % Robust equipment/tanks: partial damage
            37.5    Inf  1.00];                  %   and total loss
bands{2} = [ 5.0   12.5  0.30;                   % Sensitive equipment: partial damage
            12.5    Inf  1.00];                  %   and total loss
bands{3} = [ 5.0   12.5  0.10;                   % Buildings: light (glazing/cladding)
            12.5    Inf  1.00];                  %   and total loss (envelope ignition)
bands{4} = [37.5    Inf  1.00];                  % Fuel inventory: lost with the tank

% Trained GRNN surrogate models (m06): (x, y) -> Gamma parameters
load('./data/grnn_data.mat');

% -------------------------------------------------------------------------
% Expected Annual Loss of every asset
% -------------------------------------------------------------------------
n_assets = size(assets, 1);
E_D   = zeros(n_assets, 1);                      % Expected damage fraction [-]
L_acc = zeros(n_assets, 1);                      % Expected loss per accident [EUR]
EAL   = zeros(n_assets, 1);                      % Expected annual loss [EUR/yr]

for k = 1 : n_assets
    % Local Gamma parameters at the asset position
    a = net_a([assets{k, 2} - OFFSET, assets{k, 3} - OFFSET].');
    b = net_b([assets{k, 2} - OFFSET, assets{k, 3} - OFFSET].');

    % Expected damage fraction: sum over the bands of the stepped
    % vulnerability function of the asset class,
    %   E[D] = Sum_m d_m . [F_Gamma(I_high) - F_Gamma(I_low)]
    bnd = bands{assets{k, 5}};
    E_D(k) = 0;
    for m = 1 : size(bnd, 1)
        E_D(k) = E_D(k) + bnd(m, 3).*(gamcdf(bnd(m, 2), a, b) - ...
                                      gamcdf(bnd(m, 1), a, b));
    end

    % Expected loss per accident and Expected Annual Loss
    L_acc(k) = assets{k, 4}.*E_D(k);
    EAL(k)   = L_acc(k).*n/N;
end

% -------------------------------------------------------------------------
% Report
% -------------------------------------------------------------------------
fprintf('\n%6s %10s %10s %14s %10s %16s %14s\n', ...
        'Asset', 'x [m]', 'y [m]', 'Value [EUR]', 'E[D] [-]', ...
        'Loss/acc [EUR]', 'EAL [EUR/yr]');
for k = 1 : n_assets
    fprintf('%6s %10.1f %10.1f %14.3e %10.3e %16.3e %14.3e\n', ...
            assets{k, 1}, assets{k, 2}, assets{k, 3}, assets{k, 4}, ...
            E_D(k), L_acc(k), EAL(k));
end
fprintf('\nTotal exposed value:            %.3e EUR\n', sum([assets{:, 4}]));
fprintf('Expected loss per accident:     %.3e EUR\n', sum(L_acc));
fprintf('Total Expected Annual Loss:     %.3e EUR/yr\n', sum(EAL));

% -------------------------------------------------------------------------
% Map: assets coloured by their Expected Annual Loss (log scale),
% with the tank outlines drawn for reference
% -------------------------------------------------------------------------
fig = figure;
hold on;
% Tank outlines (circles of diameter D centred at the global positions)
theta = linspace(0, 2*pi, 100);
for t = 1 : length(D_T)
    x_circ = x_T(t) + D_T(t)/2.*cos(theta);
    y_circ = y_T(t) + D_T(t)/2.*sin(theta);
    plot(x_circ, y_circ, 'k-');
end
% Assets (marker size ~ exposed value, colour ~ EAL)
scatter([assets{:, 2}], [assets{:, 3}], ...
        300.*[assets{:, 4}]./max([assets{:, 4}]), ...
        log10(max(EAL, 1e-3)), 'filled', 'MarkerEdgeColor', 'k');
colormap('jet');
cb = colorbar;
ylabel(cb, 'log_{10} EAL [EUR/yr]');
axis equal;
grid on;
xlabel('x [m]');
ylabel('y [m]');
title('Expected Annual Loss of the assets of the complex');
saveas(fig, './figures/eal_map.fig');
saveas(fig, './figures/eal_map.png');

% -------------------------------------------------------------------------
% Save the results
% -------------------------------------------------------------------------
save('./data/eal_data.mat', 'assets', 'E_D', 'L_acc', 'EAL');
