function out = f4sb_pod_glm(varargin)
% F4SB_POD_GLM calculates POD using generalized linear model.
%
% Input: 'aData' - crack size vector [NaN] (optional)
%        'ahatData' - signal amplitude vector [NaN] (optiona)
%        'logX' - log transformation of x axis [0] (optional)
%        'logY' - log transformation of y axis [0] (optional)
%        'meanNoise' - noise mean [0.123] (optional)
%        'ahatSat' - saturation amplitude [1.505] (optional)
%        'ahatDecFactor' - factor for ahat decision [0.2] (optional)
%        'dispPlot' - plot linear model and pod [1] (optional)
% Output: out.pod - probability
%         out.aPOD50 - mean POD curve
%         out.aPOD95 - 95% confidence bound of POD curve
%         out.a50 - crack lengths of 50% POD
%         out.a90 - crack lengths of 90% POD
%         out.a9095 - crack lengths of 90% POD with 95% confidence bound
%
% Revision history
% 061016 LDY Code was modified from 'tang_fcg_pod_glm_v2.m' for J85 comp.
%            1st stage rotor blade POD calculation.

% Set default.
aData = NaN;           % Crack size
ahatData = NaN;        % Signal amplitude
logX = 0;              % Log transformation of x axis
logY = 0;              % Log transformation of y axis
ampNoise = 0.123;      % Noise amplitude, [dB]
ampSat = 1.505;        % Saturation amplitude, [dB]
ampDecFactor = 0.2;    % Factor for decision amplitude
dispPlot = 1;          % Plot linear model and pod.

% Input arguments
args = varargin;
for i = 1:2:length(args)
    switch args{i}
        case 'aData', aData = args{i+1};
        case 'ahatData', ahatData = args{i+1};
        case 'logX', logX = args{i+1};
        case 'logY', logY = args{i+1};
        case 'meanNoise', ampNoise = args{i+1};
        case 'ahatSat', ampSat = args{i+1};
        case 'ahatDecFactor', ampDecFactor = args{i+1};
        case 'dispPlot', dispPlot = args{i+1};
        otherwise, error(['invalid argument name: ' args{i}]);
    end
end

% Constant
X = 0:0.1:4;    % For plot
if logX
    X = log(X);
end

% Get decision amplitude, [dB]
ahatDec = (ampSat - ampNoise)*ampDecFactor + ampNoise;
if logY
    ahatDec = log(ahatDec);
end

% Data --------------------------------------------------------------------
% Read data.
if isnan(aData)
    data = xlsread('a_vs_ahat', 1);
    aData = data(:, 1);
    ahatData = data(:, 2);
end

% Indices for noise, signal, saturation
idxNoise = ahatData <= ampNoise;
idxSat = ahatData >= ampSat;
idxSignal = logical(~idxNoise.*~idxSat);

% ##### Noise #####
% Crack length
aNoise = aData(idxNoise);
if logX
    aNoise = log(aNoise);
end

% Amplitude
ahatNoise = ahatData(idxNoise);
if logY
    % Substitute zeros to prevent -inf which causes failure in glm creation.
    idxZeroNoise = ahatNoise==0;
    ahatNoise(idxZeroNoise) = eps;
    
    ahatNoise = log(ahatNoise);
    ampNoise = log(ampNoise);
end

% ##### Signal #####
nSignalData = sum(idxSignal);

% Crack length
aSignal = aData(idxSignal);
if logX
    aSignal = log(aSignal);
end

% Amplitude
ahatSignal = ahatData(idxSignal);
if logY
    % Substitute zeros to prevent -inf which causes failure in glm creation.
    idxZeroSignal = ahatSignal==0;
    ahatSignal(idxZeroSignal) = eps;
    
    ahatSignal = log(ahatSignal);
end

% Get averages of a for signal.
aSignalBar = mean(aSignal);

% ##### Saturation #####
% Crack length
aSat = aData(idxSat);
if logX
    aSat = log(aSat);
end

% Amplitude
ahatSat = ahatData(idxSat);
if logY
    ahatSat = log(ahatSat);
    ampSat = log(ampSat);
end

% Get linear model. -------------------------------------------------------
[estimator, logRMSE, CovTobit] = TOBIT( ...
    [ahatSignal; ahatNoise; ahatSat], ...
    [aSignal; aNoise; aSat], ...
    'lb', ampNoise, ...    % Lower bound to censoring
    'ub', ampSat, ...       % Upper bound to censoring
    'add_constant', 1, ...
    'verbose', dispPlot);

b0 = estimator(1);      % Intercept
b1 = estimator(2);      % Slope
RMSE = exp(logRMSE);    % Root mean square error of residuals

% Degrees of freedom
df = nSignalData - 2;

% Sum of squares
SSx = sum((aSignal - aSignalBar).^2);

% Linear model confidence bounds
tCrit = abs(tinv(0.025, df));

% Confidence interval for mu_Y|X (the mean of the population of Y values
% corresponding to Xi)
seMean = RMSE*sqrt(1/nSignalData + (X - aSignalBar).^2/SSx);
lbMean = (b0 + b1*X) - tCrit*seMean;    % 95% lower bound
ubMean = (b0 + b1*X) + tCrit*seMean;    % 95% lower bound

% Confidence interval for Yi (an individual predicted value)
seYi = RMSE*sqrt(1 + 1/nSignalData + (X - aSignalBar).^2/SSx);
lbYi =(b0 + b1*X) - tCrit*seYi;    % 95% lower bound
ubYi =(b0 + b1*X) + tCrit*seYi;    % 95% upper bound

% POD model ---------------------------------------------------------------
% Get POD model parameters.
muGLM =(ahatDec - b0)/b1;
sigmaGLM = RMSE/b1;

% Covariance matrix
phi =-1/b1*[1, 0; muGLM, sigmaGLM; 0, -RMSE];
CovPOD = phi'*CovTobit*phi;

% Mean POD curve
pod = 0:0.001:1;
aPOD50 = muGLM + norminv(pod)*sigmaGLM;

% Standard error for confidence bounds
sePOD = sqrt(CovPOD(1) + norminv(pod).^2*CovPOD(4) + ...
    2*norminv(pod)*CovPOD(2));

% 95% Confidence bound
aPOD95 = aPOD50 + norminv(0.95)*sePOD;

% Transform if necessary.
if logX
    aPOD50 = exp(aPOD50);
    aPOD95 = exp(aPOD95);
end

% Get meaningful crack sizes.
a5050 = muGLM + norminv(0.5)*sigmaGLM;
a9050 = muGLM + norminv(0.9)*sigmaGLM;
a9095 = a9050 + norminv(0.95)*sePOD(pod==0.9);
if logX
    a5050 = exp(a5050);
    a9050 = exp(a9050);
    a9095 = exp(a9095);
end

% Return ------------------------------------------------------------------
out.pod = pod;
out.aPOD50 = aPOD50;    % Mean POD curve
out.aPOD95 = aPOD95;    % 95% confidence bound
out.a5050 = a5050;
out.a9050 = a9050;
out.a9095 = a9095;

% Plot
if dispPlot
    % Plot 1 --------------------------------------------------------------
    figure(1)
    if logX && logY
        h(1) = loglog(exp(X), exp(b0+b1*X), '-r', 'LineWidth', 2); hold on
        h(2) = loglog(exp(X), exp(lbMean), '--b', 'LineWidth', 1);
        h(3) = loglog(exp(X), exp(ubMean), '--b', 'LineWidth', 1);
        h(4) = loglog(exp(X), exp(lbYi), ':b', 'LineWidth', 1);
        h(5) = loglog(exp(X), exp(ubYi), ':b', 'LineWidth', 1);
        h(6) = loglog(exp(aSignal(~idxZeroSignal)), exp(ahatSignal(~idxZeroSignal)), '.');
        h(7) = loglog(exp(aNoise(~idxZeroNoise)), exp(ahatNoise(~idxZeroNoise)), '.r');
        h(8) = loglog(exp(aSat), exp(ahatSat), '.r');
        line(exp(X), exp(ampNoise)*(ones(size(X))), 'Color', 'k', 'LineStyle', '--');
        line(exp(X), exp(ahatDec)*ones(size(X)), 'Color', 'k', 'LineStyle', '--');
        line(exp(X), exp(ampSat)*ones(size(X)), 'Color', 'k', 'LineStyle', '--'); hold off
        
        axisRange = [0, inf, 0, 2];
    else
        h(1) = plot(X, b0+b1*X, '-r', 'LineWidth', 2); hold on
        h(2) = plot(X, lbMean, '--b', 'LineWidth', 1);
        h(3) = plot(X, ubMean, '--b', 'LineWidth', 1);
        h(4) = plot(X, lbYi, ':b', 'LineWidth', 1);
        h(5) = plot(X, ubYi, ':b', 'LineWidth', 1);
        h(6) = plot(aSignal, ahatSignal, '.');
        h(7) = plot(aNoise, ahatNoise, '.r');
        h(8) = plot(aSat, ahatSat, '.r');
        line(X, ampNoise(ones(size(X))), 'Color', 'k', 'LineStyle', '--');
        line(X, ahatDec*ones(size(X)), 'Color', 'k', 'LineStyle', '--');
        line(X, ampSat*ones(size(X)), 'Color', 'k', 'LineStyle', '--'); hold off
        
        axisRange = [0, inf, 0, 2];
    end
    
    % Format
    grid on
    axis(axisRange);
    set(gca, 'FontSize', 14)
    set(gca, 'FontWeight', 'bold')
    
    % legend
    %h = get(gca, 'Children');
    leg = legend(h(1:2), {'Linear Model', 'Confidence Bound'});
    set(leg, 'Location', 'NorthWest');
    set(leg, 'Color', 'none');
    set(leg, 'FontSize', 14);
    set(leg, 'FontWeight', 'bold');
    
    % Label
    xl = xlabel('Crack Length (mm)');
    set(xl, 'FontSize', 15);
    set(xl, 'FontWeight', 'bold');
    yl = ylabel('Amplitude (dB)');
    set(yl, 'FontSize', 15);
    set(yl, 'FontWeight', 'bold');
    
    % Plot 2 --------------------------------------------------------------
    figure(2)
    plot(aPOD50, pod, '-r', 'LineWidth', 2); hold on
    plot(aPOD95, pod, '--', 'LineWidth', 1); hold off
    text(aPOD95(end-1)*0.8, 0.70,['a_{50/50} = ', num2str(a5050, '%2.3f')]);
    text(aPOD95(end-1)*0.8, 0.62,['a_{90/50} = ', num2str(a9050, '%2.3f')]);
    text(aPOD95(end-1)*0.8, 0.54,['a_{90/95} = ', num2str(a9095, '%2.3f')]);
    text(aPOD95(end-1)*0.8, 0.46,['$\hat{\mu}$ = ', num2str(muGLM, '%2.3f')], 'Interpreter', 'Latex');
    text(aPOD95(end-1)*0.8, 0.38,['$\hat{\sigma}$ = ', num2str(sigmaGLM, '%2.3f')], 'Interpreter', 'Latex');
    
    % Format
    grid on
    set(gca, 'FontSize', 14);
    set(gca, 'FontWeight', 'bold');
    
    % legend
    h = get(gca, 'Children');
    leg = legend([h(7), h(6)], {'Mean POD', '95% Lower Bound'});
    set(leg, 'Location', 'NorthWest');
    set(leg, 'Color', 'none');
    set(leg, 'FontSize', 14);
    set(leg, 'FontWeight', 'bold');
    
    % Label
    xl = xlabel('Crack Length (mm)');
    set(xl, 'FontSize', 15);
    set(xl, 'FontWeight', 'bold');
    yl = ylabel('Probability of Detection');
    set(yl, 'FontSize', 15);
    set(yl, 'FontWeight', 'bold');
end

end