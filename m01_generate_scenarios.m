% =========================================================================
% Paper: Towards a Digital Twin for the Quantitative Risk Assessment of
%     Domino Effects in Fuel Storage Facilities
% Repository: https://github.com/ramon-quiza/qra-domino-fuel-storage
% File: Generates, randomly, the accident scenarios
% -------------------------------------------------------------------------
% (c) 2026, Centre for Advanced and Sustainable Manufacturing Studies
%           University of Matanzas
% =========================================================================
clear all;                                       % Clear workspace variables
clc;                                             % Clear command window

% -------------------------------------------------------------------------
% Simulation parameters
% -------------------------------------------------------------------------
% Accident likelihood: probability that a given tank suffers a primary
% accident in a scenario. The value of 2e-5 corresponds to the accident
% frequency of fixed roof atmospheric storage tanks (~2e-5 per tank per
% year) reported in:
%   Taveau, J. (2011). Explosion of fixed roof atmospheric storage
%   tanks, part 1: Background and review of case histories. Process
%   Safety Progress, 30(4), 381-392. https://doi.org/10.1002/prs.10459
% Note: since this frequency is expressed per tank-year, each Monte
% Carlo trial represents one year of operation of the facility. The
% value is interpreted as the frequency of a full-surface tank fire with
% immediate ignition (conservative envelope, see Section 2.1 of the
% paper). Every annual indicator of the pipeline scales linearly with it.
ACCIDENT_LH = 2e-5;
SCENARIOS_COUNT = 1e8;                           % Number of Monte Carlo trials (facility-years)
% Note: no random seed is fixed, so successive runs give statistically
% equivalent but not identical samples. Add rng(<seed>) here to
% reproduce a given run exactly.

% -------------------------------------------------------------------------
% Tank farm layout: position and geometry of each storage tank
% (one column per tank; 5 tanks in total), in the GLOBAL coordinate
% system of the complex (see the layout figure and data tables)
% -------------------------------------------------------------------------
x =   [150.0  220.0  177.0  162.5  217.5];       % X-axis centre position [m]
y =   [150.0  150.0  177.0  240.0  240.0];       % Y-axis centre position [m]
D =   [ 15.0   10.0   10.0   22.0   22.0];       % Diameter [m]
H =   [ 12.0    8.0   10.0   12.0   12.0];       % Height [m]

% -------------------------------------------------------------------------
% Thermophysical properties of the fuel contained in the tanks
% - Tank 1: Light naphtha (n-hexane) [https://www.chemeo.com/cid/14-924-0/n-Hexane]
% - Tanks 2 & 3: Jet fuel (n-decane) [https://www.chemeo.com/cid/44-644-8/Decane]
% - Tanks 4 & 5: Diesel (n-dodecane) [https://www.chemeo.com/cid/34-125-5/Dodecane]
% -------------------------------------------------------------------------
T_b = [ 68.7  174.2  174.2  216.3  216.3];       % Boiling point [oC]
H_c = [ 44.7   47.6   47.6   47.5   47.5].*1e3;  % Heat of combustion [kJ/kg]
H_v = [366.0  361.0  361.0  256.0  256.0];       % Heat of vaporization [kJ/kg]
c_P = [  2.27   2.21   2.21   2.21   2.21];      % Specific heat [kJ/(kg.K)]

% -------------------------------------------------------------------------
% Monte Carlo sampling of primary (initiating) accidents
% -------------------------------------------------------------------------
% Uniform random numbers in [0, 1): one row per scenario, one column per tank
r = rand(SCENARIOS_COUNT, length(D));

% Row indices of the scenarios where at least one tank suffers a primary
% accident (r < ACCIDENT_LH); all other scenarios are discarded
id = find(sum(r < ACCIDENT_LH, 2) > 0);

% Keep only the relevant scenarios and encode them as a 0/1 matrix:
%   id(i, j) = 1 -> tank j suffers a primary accident in scenario i
%   id(i, j) = 0 -> tank j is initially intact in scenario i
% Later, id(i, j) will also store the escalation step (c) at which an
% initially intact tank is ignited by thermal radiation (domino effect)
id = double(r(id, :) < ACCIDENT_LH);

% -------------------------------------------------------------------------
% Domino-effect (escalation) simulation for each accident scenario
% -------------------------------------------------------------------------
for i = 1 : length(id)                           % Loop over accident scenarios
    fl_stop = true;                              % Flag: keep iterating while new escalations occur
    c = 1;                                       % Escalation step counter (1 = primary accident)
    while fl_stop
        % disp(c);
        c = c + 1;                               % Advance to the next escalation step
        % Heat released by each burning tank (pool fire model);
        % intact tanks (id = 0) release no heat
        for j = 1 : length(D)
            Q_R(i, j) = min(id(i, j), 1.0).*poolfire(D(j), T_b(j), H_c(j), H_v(j), c_P(j));
        end
        % Thermal radiation intensity received by each tank from all
        % the burning tanks in the facility
        for j = 1 : length(D)
            I_R(i, j) = radiation(j, x, y, D, Q_R(i, :));
        end
        % Escalation trial: each intact tank (Q_R <= 0) fails - and is
        % assumed to ignite - with the probit probability of the
        % radiation it currently receives (Bernoulli trial, see
        % escalation.m); trials at successive steps are independent
        fl_stop = false;
        for j = 1 : length(D)
            if ((Q_R(i, j) <= 0.0) && escalation(I_R(i, j), D(j), H(j)))
                id(i, j) = c;                    % Record the step at which tank j ignites
                %disp([c, id(i, :)]);
                fl_stop = true;                  % A new fire appeared: iterate again
            end
        end
        
    end
end

% -------------------------------------------------------------------------
% Save the results:
%   id  - ignition step of each tank (0 = intact, 1 = primary, >1 = domino)
%   Q_R - heat released by each tank in each scenario [kW]
%   I_R - radiation intensity received by each tank in each scenario [kW/m2]
% The number of retained scenarios, n = size(id, 1), is the value that
% m07-m10 need as the numerator of the annual accident probability n/N.
% -------------------------------------------------------------------------
save('./data/scenarios.mat', 'id', 'Q_R', 'I_R');