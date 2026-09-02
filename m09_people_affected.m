% =========================================================================
% Paper: Towards a Digital Twin for the Quantitative Risk Assessment of
%     Domino Effects in Fuel Storage Facilities
% Repository: https://github.com/ramon-quiza/qra-domino-fuel-storage
% File: Evaluates the expected number of people affected (fatalities)
%     in the buildings surrounding the facility
% -------------------------------------------------------------------------
% (c) 2026, Centre for Advanced and Sustainable Manufacturing Studies
%           University of Matanzas
% =========================================================================
% Given the buildings of the complex (see the data tables), each with
% its location and its usual number of users, this script uses the
% digital twin (GRNN surrogates of the local Gamma distribution of the
% radiation intensity, m06) and the vulnerability models implemented in
% people_affected to estimate the expected number of people affected in
% each building, both conditional on the occurrence of an accident and
% on an annual basis (multiplying by the annual accident probability
% n/N, see m07). A fraction f_out of the users is assumed to be
% outdoors, and the remainder sheltered inside the buildings (see
% people_affected for the vulnerability models and their references).
% All positions are given in the GLOBAL coordinate system of the
% complex, consistent with the rest of the pipeline (m01-m08).
% =========================================================================
clear all;                                       % Clear workspace variables
clc;                                             % Clear command window

% -------------------------------------------------------------------------
% Tank farm layout (GLOBAL frame; needed only to draw the map)
% -------------------------------------------------------------------------
x =   [150.0  220.0  177.0  162.5  217.5];       % X-axis centre position [m]
y =   [150.0  150.0  177.0  240.0  240.0];       % Y-axis centre position [m]
D =   [ 15.0   10.0   10.0   22.0   22.0];       % Diameter [m]

% -------------------------------------------------------------------------
% Annual accident probability (see m01 and m07)
% -------------------------------------------------------------------------
N = 1e8;                                         % Monte Carlo trials generated in m01 (one per facility-year)
n = 9908;                                        % Trials that resulted in an accident: MUST equal size(id, 1) of the
                                                 % scenarios.mat actually used (n = 10048 in the run reported in the paper)

% -------------------------------------------------------------------------
% Population model
% -------------------------------------------------------------------------
f_out = 0.05;                                    % Fraction of the users outdoors [-] (time-averaged value; Bevi prescribes 0.07 by day and 0.01 by night)
t_exp = 20;                                      % Effective exposure time of outdoor people [s] (20-30 s recommended, see people_affected)

% -------------------------------------------------------------------------
% Buildings of the complex (GLOBAL coordinates) and their usual users:
% B1 administrative office building; B2 heavy engineering works;
% B3-B8 2-storey detached houses; B9-B12 warehouses;
% B13-B22 5-storey blocks of flats (see the data tables)
% -------------------------------------------------------------------------
%        B1     B2     B3     B4     B5     B6     B7     B8     B9    B10    B11    B12    B13    B14    B15    B16    B17    B18    B19    B20    B21    B22
x_B = [215.3  159.1  292.0  317.0  342.0  367.0   96.9   71.9  327.5  337.5  337.5   83.7  338.8  338.8  338.8  248.8  248.8  248.8  158.8  158.8  158.8   66.0];
y_B = [304.0  374.2  369.5  369.5  369.5  369.5  369.5  369.5  276.6  204.6  159.6  250.0   80.6   45.6   10.6   80.6   45.6   10.6   80.6   45.6   10.6   21.5];
P_B = [   23     44      4      4      4      4      4      4      8      8      8      8     80     80     80     80     80     80     80     80     80     80];   % Usual (24-h averaged) occupancy

% Trained GRNN surrogate models (m06): (x, y) -> Gamma parameters
load('./data/grnn_data.mat');

% -------------------------------------------------------------------------
% Expected number of people affected in each building
% -------------------------------------------------------------------------
% Conditional on the occurrence of an accident
[N_acc, p_out, p_in] = people_affected(x_B, y_B, P_B, ...
                                       net_a, net_b, f_out, t_exp);
% On an annual basis (accident probability n/N per facility-year)
N_yr = N_acc.*n/N;

% -------------------------------------------------------------------------
% Report
% -------------------------------------------------------------------------
fprintf('\n%8s %10s %10s %10s %12s %14s %16s\n', ...
        'Building', 'x [m]', 'y [m]', 'Users', ...
        'p_out [-]', 'p_in [-]', 'Affected/acc');
for k = 1 : length(P_B)
    fprintf('%8s %10.1f %10.1f %10d %12.3e %14.3e %16.3e\n', ...
            sprintf('B%d', k), x_B(k), y_B(k), P_B(k), ...
            p_out(k), p_in(k), N_acc(k));
end
fprintf('\nTotal users of the complex:            %d\n', sum(P_B));
fprintf('Expected people affected per accident: %.3e\n', sum(N_acc));
fprintf('Expected people affected per year:     %.3e\n', sum(N_yr));

% -------------------------------------------------------------------------
% Map: facility layout and buildings (marker size ~ users,
% colour ~ log10 of the expected number of people affected per year)
% -------------------------------------------------------------------------
fig = figure;
hold on;
% Tank outlines (circles of diameter D centred at (x, y))
theta = linspace(0, 2*pi, 100);
for t = 1 : length(D)
    x_circ = x(t) + D(t)/2.*cos(theta);
    y_circ = y(t) + D(t)/2.*sin(theta);
    plot(x_circ, y_circ, 'k-');
end
% Buildings
scatter(x_B, y_B, 10 + 5.*P_B, log10(N_yr), 'filled', 'MarkerEdgeColor', 'k');
colormap('jet');
cb = colorbar;
ylabel(cb, 'log_{10} expected people affected per year');
axis equal;
grid on;
xlabel('x [m]');
ylabel('y [m]');
set(gca, 'XTick', [-50 : 50 : 450]);
set(gca, 'XLim', [-50, 450]);
set(gca, 'YTick', [-50 : 50 : 450]);
set(gca, 'YLim', [-50, 450]);
title(sprintf('People affected (f_{out} = %.0f %%, t_{exp} = %d s)', ...
              100*f_out, t_exp));
saveas(fig, './figures/people_affected.fig');
saveas(fig, './figures/people_affected.png');

% -------------------------------------------------------------------------
% Save the results
% -------------------------------------------------------------------------
save('./data/people_data.mat', 'x_B', 'y_B', 'P_B', ...
     'p_out', 'p_in', 'N_acc', 'N_yr');
