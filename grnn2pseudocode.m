% =========================================================================
% Paper: Towards a Digital Twin for the Quantitative Risk Assessment of
%     Domino Effects in Fuel Storage Facilities
% Repository: https://github.com/ramon-quiza/qra-domino-fuel-storage
% File: Writes the evaluation pseudocode of a trained GRNN, including
%     the fitted parameter values, so that the surrogate model can be
%     documented and re-implemented on any platform
% -------------------------------------------------------------------------
% (c) 2026, Centre for Advanced and Sustainable Manufacturing Studies
%           University of Matanzas
% =========================================================================
% Inputs:
%   net        - Trained GRNN, as returned by newgrnn (m06)
%   file_name  - Path of the text file where the pseudocode is written
% Output:
%   err_max    - Maximum absolute difference between the network output
%                and the transcribed pseudocode formula, evaluated on a
%                set of random verification points (a value close to
%                machine precision certifies that the pseudocode is a
%                faithful transcription of the trained network)
% -------------------------------------------------------------------------
% Model transcribed: a GRNN [1] is a normalised radial-basis-function
% estimator. MATLAB's newgrnn implementation stores:
%   - the input weights IW{1,1}: the training inputs (one centre c_i
%     per neuron/training sample);
%   - the layer weights LW{2,1}: the training targets w_i;
%   - the bias b{1}: identical for every neuron and equal to
%     0.8326/spread, so that a neuron responds 0.5 when the input lies
%     at a distance equal to the spread from its centre.
% The network output is
%   y(x) = sum_i w_i.exp(-(b.d_i)^2) / sum_i exp(-(b.d_i)^2),
%   d_i = ||x - c_i||,
% which is equivalent to the Nadaraya-Watson form of the paper,
% y = sum w_i.exp(-d_i^2/(2.sigma^2)) / sum exp(-d_i^2/(2.sigma^2)),
% with sigma = spread/1.1774.
% -------------------------------------------------------------------------
% References:
%   [1] Specht, D.F. (1991). A general regression neural network. IEEE
%       Transactions on Neural Networks, 2(6), 568-576.
%       https://doi.org/10.1109/72.97934
% =========================================================================
function err_max = grnn2pseudocode(net, file_name)

    % ---------------------------------------------------------------------
    % Extraction of the fitted parameters from the network object
    % ---------------------------------------------------------------------
    C = net.IW{1, 1};                            % Centres [Q x R]: training inputs
    w = net.LW{2, 1};                            % Weights [1 x Q]: training targets
    b = net.b{1}(1);                             % Bias [-] (identical for all neurons)
    Q = size(C, 1);                              % Number of neurons (training samples)
    R = size(C, 2);                              % Input dimension
    spread = 0.8326./b;                          % Spread of the radial basis functions
    sigma  = spread./1.1774;                     % Equivalent Gaussian sigma (paper notation)

    % ---------------------------------------------------------------------
    % Pseudocode writing
    % ---------------------------------------------------------------------
    fid = fopen(file_name, 'w');
    fprintf(fid, '=====================================================================\n');
    fprintf(fid, ' Evaluation pseudocode of the trained GRNN surrogate model\n');
    fprintf(fid, ' Generated automatically by grnn2pseudocode on %s\n', datestr(now));
    fprintf(fid, '=====================================================================\n\n');

    fprintf(fid, 'CONSTANTS:\n');
    fprintf(fid, '    Q      = %d          // number of neurons (training samples)\n', Q);
    fprintf(fid, '    R      = %d          // input dimension\n', R);
    fprintf(fid, '    spread = %.10g       // spread of the radial basis functions\n', spread);
    fprintf(fid, '    b      = 0.8326/spread = %.10g\n', b);
    fprintf(fid, '    // equivalent Gaussian sigma (Nadaraya-Watson form): %.10g\n\n', sigma);

    fprintf(fid, 'ALGORITHM Evaluate_GRNN(x[1..R]) RETURNS y:\n');
    fprintf(fid, '    num <- 0\n');
    fprintf(fid, '    den <- 0\n');
    fprintf(fid, '    FOR i <- 1 TO Q DO\n');
    fprintf(fid, '        // Euclidean distance to the centre of neuron i\n');
    fprintf(fid, '        d <- sqrt( sum_{j=1..R} ( x[j] - c[i][j] )^2 )\n');
    fprintf(fid, '        // Radial basis response of neuron i\n');
    fprintf(fid, '        a <- exp( -( b * d )^2 )\n');
    fprintf(fid, '        num <- num + w[i] * a\n');
    fprintf(fid, '        den <- den + a\n');
    fprintf(fid, '    END FOR\n');
    fprintf(fid, '    y <- num / den\n');
    fprintf(fid, 'END ALGORITHM\n\n');

    fprintf(fid, 'FITTED PARAMETERS (one row per neuron):\n');
    % Header of the parameter table
    fprintf(fid, '%8s', 'i');
    for j = 1 : R
        fprintf(fid, '%16s', sprintf('c[i][%d]', j));
    end
    fprintf(fid, '%16s\n', 'w[i]');
    % Parameter values (centres and weights)
    for i = 1 : Q
        fprintf(fid, '%8d', i);
        fprintf(fid, '%16.8g', C(i, :));
        fprintf(fid, '%16.8g\n', w(i));
    end

    % ---------------------------------------------------------------------
    % Verification: the transcribed formula is evaluated on random
    % points of the training domain and compared against the network
    % ---------------------------------------------------------------------
    N_CHECK = 20;                                % Number of verification points
    x_chk = min(C) + rand(N_CHECK, R).*(max(C) - min(C));
    err_max = 0;
    for k = 1 : N_CHECK
        d2 = sum((C - x_chk(k, :)).^2, 2);       % Squared distances to the centres
        a  = exp(-(b.^2).*d2);                   % Radial basis responses
        y_pseudo = (w*a)./sum(a);                % Pseudocode formula
        y_net    = net(x_chk(k, :).');           % Network evaluation
        err_max  = max(err_max, abs(y_pseudo - y_net));
    end
    fprintf(fid, '\nVERIFICATION: max |y_pseudocode - y_network| = %.3e over %d random points\n', ...
            err_max, N_CHECK);
    fclose(fid);

    fprintf('Pseudocode written to %s (Q = %d neurons; max verification error = %.3e)\n', ...
            file_name, Q, err_max);
end
