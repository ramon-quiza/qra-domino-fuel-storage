% =========================================================================
% Paper: Towards a Digital Twin for the Quantitative Risk Assessment of
%     Domino Effects in Fuel Storage Facilities
% Repository: https://github.com/ramon-quiza/qra-domino-fuel-storage
% File: Consolidate radiation data for further analysis
% -------------------------------------------------------------------------
% (c) 2026, Centre for Advanced and Sustainable Manufacturing Studies
%           University of Matanzas
% =========================================================================
% This script merges the radiation fields of the 32 fire configurations
% computed by m02 into a single, flat dataset of (x, y, I) samples. Each
% configuration is replicated as many times as it occurred in the Monte
% Carlo simulation (C), so that the frequency of each configuration is
% embodied in the dataset itself. The rows are finally shuffled at
% random, leaving the dataset ready for further analysis (e.g. the
% training of a data-driven model of the radiation distribution).
% Note: this script is auxiliary and is not required by the main
% pipeline (m01 -> m02 -> m04 -> m05 -> m06 -> m07-m10); the variable
% "data" is left in the workspace and is not saved to disk.
% =========================================================================
clear all;                                       % Clear workspace variables
clc;                                             % Clear command window

% -------------------------------------------------------------------------
% Consolidation of the per-configuration radiation fields
% -------------------------------------------------------------------------
% Accumulators for the flattened samples: point coordinates [m] and
% received radiation intensity [kW/m2]
X_data = []; Y_data = []; Z_data = [];

% Data files produced by m02 (one per fire configuration), each holding:
%   X, Y - grid coordinates [m]; Z - radiation field [kW/m2]
%   C    - occurrence count;     S - configuration label
data_files = dir('./data/scenarios/*.mat');

C_total = 0;                                     % Total number of consolidated scenarios
for i = 1 : length(data_files)
    load(['./data/scenarios/', data_files(i).name])
    % Only configurations observed in the Monte Carlo simulation
    % contribute to the dataset; note that the empty configuration
    % ('00000') never contributes, since m01 retains only scenarios
    % with at least one burning tank
    if C > 0
        % Flatten the (N+1) x (N+1) grids into column vectors and
        % replicate them C times, so that each configuration is
        % weighted by its frequency of occurrence
        X_data = [X_data; repmat(reshape(X, numel(X), 1), C, 1)];
        Y_data = [Y_data; repmat(reshape(Y, numel(Y), 1), C, 1)];
        Z_data = [Z_data; repmat(reshape(Z, numel(Z), 1), C, 1)];

        C_total = C_total + C;
    end
end

% -------------------------------------------------------------------------
% Assembly and random shuffling of the dataset
% -------------------------------------------------------------------------
% One sample per row: [x-coordinate, y-coordinate, radiation intensity]
data = [X_data, Y_data, Z_data];

% Shuffle the rows at random to remove the ordering by configuration
% (advisable before splitting the dataset for training/validation)
rnd_rows = randperm(size(data, 1));
data = data(rnd_rows, :);

% Quick verification plot of the last loaded configuration (the
% intensity is scaled by the grid size only for display purposes)
fig = figure;
plot3(X, Y, Z./length(Z), 'b.');
saveas(fig, './figures/plot3d.fig');
saveas(fig, './figures/plot3d.png');