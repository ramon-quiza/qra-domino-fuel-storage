% =========================================================================
% Paper: Towards a Digital Twin for the Quantitative Risk Assessment of
%     Domino Effects in Fuel Storage Facilities
% Repository: https://github.com/ramon-quiza/qra-domino-fuel-storage
% File: Organizes the radiation fields into a scenario-wise 3D array
%     for the distribution fitting
% -------------------------------------------------------------------------
% (c) 2026, Centre for Advanced and Sustainable Manufacturing Studies
%           University of Matanzas
% =========================================================================
% This script rearranges the radiation fields of the fire configurations
% computed by m02 into a single 3D array Z_data(n, i, j), where n indexes
% the scenario occurrence and (i, j) the grid point. Each configuration
% is replicated as many times as it occurred in the Monte Carlo
% simulation (C), so that, at every grid point, the first dimension of
% Z_data holds a frequency-weighted sample of the radiation intensity,
% ready for the point-wise distribution fitting performed in m05.
% =========================================================================
clear all                                        % Clear workspace variables
clc                                              % Clear command window

% -------------------------------------------------------------------------
% Consolidation of the per-configuration radiation fields
% -------------------------------------------------------------------------
% Data files produced by m02 (one per fire configuration), each holding:
%   X, Y - grid coordinates [m]; Z - radiation field [kW/m2]
%   C    - occurrence count;     S - configuration label
data_files = dir('./data/scenarios/*.mat');

n = 1;                                           % Scenario occurrence counter
for i = 1 : length(data_files)
    load(['./data/scenarios/', data_files(i).name])
    % Replicate the radiation field of this configuration C times, so
    % that its frequency of occurrence is embodied in the sample; the
    % empty configuration ('00000') never contributes, since m01
    % retains only scenarios with at least one burning tank (C = 0)
    for j = 1 : C
        Z_data(n, :, :) = Z;
        n = n + 1;
    end
    % Progress report
    fprintf('Organised %d of %d\n', i, length(data_files));
end

% -------------------------------------------------------------------------
% Save the scenario-wise radiation sample:
%   Z_data - radiation intensity [kW/m2], dimensions:
%            (scenario occurrence) x (grid row) x (grid column)
% -------------------------------------------------------------------------
save('./data/z_data.mat', 'Z_data');
