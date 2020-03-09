function [estimator, sigma, cov_Hessian, ME1, ME2]= TOBIT(Y, X, varargin)
% TOBIT estimates Tobit model (censored dependent variable) and marginal
% effects.
% ----------------------------------
% Model:
% Yi* = Xi * Beta + ui , where ui ~ N(0,s^2)
% Yi* is unobservable. We observe sensored version Yi with lb and ub
% ----------------------------------
% Algorithm:
% Maximum likelihood,
% Censored part follows normal c.d.f., uncensored part has normal p.d.f.
% ----------------------------------
% Usage:
% Y = dependent variable (n * 1 vector)
% X = regressors (n * k matrix)
% 'lb' = lower bound due to censoring, support inf, NaN [NaN] (optional)
% 'ub' = upper bound due to censoring, support inf, NaN [NaN] (optional)
% 'add_constant' = whether to add a constant to X [1] (optional)
% 'verbose' = display censoring bounds and log likelihood [1] (optional)
% ----------------------------------
% Returns:
% estimator = estimator corresponding to the k regressors
% sigma = estimated standard deviation of disturbances
% cov_Hessian = covariance matrix of the estimator
% ME1 = marginal effects (average data)
% ME2 = marginal effects (individual average)
%
% Written by Hang Qian, Iowa State University
% Contact me: matlabist@gmail.com 
%
% Revision history
% 060816 LDY Tobit model was parameterized in terms of the log of the
%            variance.
% 030620 LDY 'verbose' was added to suppress printing unecessary
%            information.

% Set default.
lb = NaN;            % Lower bound due to censoring
ub = NaN;            % Upper bound due to censoring
add_constant = 0;    % Whether to add a constant to X 
verbose = 1;         % Display censoring bounds and log likelihood

% Input arguments
args = varargin;
for i = 1:2:length(args)
    switch args{i}
        case 'lb', lb = args{i+1};
        case 'ub', ub = args{i+1};
        case 'add_constant', add_constant = args{i+1};
        case 'verbose', verbose = args{i+1};
        otherwise, error(['invalid argument name: ' args{i}]);
    end
end

% Check the installation of statistics toolbox.
try
    JustTest = normcdf(1);
catch
    disp('')
    disp('Oooopse, Matlab statistics toolbox is not installed.')
    disp('You may download a compatibility package on my website.')
    disp('http://www.public.iastate.edu/~hqi/toolkit')
    error('Program exits.')
end

[nrow_x,ncol_x] = size(X);
[nrow_y,ncol_y] = size(Y);
if nrow_x < ncol_x; X = X'; ncol_x = nrow_x; end
if nrow_y < ncol_y; Y = Y'; ncol_y = nrow_y; end
if ncol_x < ncol_y; Y_temp = Y; Y = X; X = Y_temp; end

[nobs, nreg]= size(X);
if add_constant == 1
    if verbose
        disp('A constant is added to X')
    end
    X = [ones(nobs,1),X];
    nreg = nreg + 1;
end

% Set lower and upper bounds
if isnan(lb) && isnan(ub)
    if verbose
        disp('Censor points are not specified, the program will assign one according to thedata.')
    end
    
    lb = min(Y);
    if sum(Y==lb) < nobs*0.01
        lb = -inf;
    end
    
    ub = max(Y);
    if sum(Y==ub) < nobs*0.01
        ub = inf;
    end
    
    if verbose
        disp(['The censored lower bound is assigned asY = ',num2str(lb)])
        disp(['The censored upper bound is assigned as: Y = ',num2str(ub)])
    end
elseif ~isnan(lb) && isnan(ub)
    ub = inf;
    if verbose
        disp(['User specified censored lower bound ',num2str(lb),', without upper bound'])
    end
elseif isnan(lb) && ~isnan(ub)
    lb = -inf;
    if verbose
        disp(['User specified censored lower bound ', num2str(lb)])
    end
else
    if verbose
        disp(['User specified censored lower bound ',num2str(lb),', and upper bound ',num2str(ub)])
    end
end

% Get indicies for censored elements.
lb_index = (Y <= lb);
X_lb = X(lb_index,:);

ub_index = (Y >= ub);
X_ub = X(ub_index,:);

continuous_index = ~(lb_index | ub_index);
Y_continous = Y(continuous_index);
X_continous = X(continuous_index,:);
num_continous = length(Y_continous);

coeff_OLS = (X'*X)\(X'*Y);
sigma_OLS = sqrt((Y-X*coeff_OLS)'*(Y-X*coeff_OLS)/nobs);
retain_portion = normcdf(ub,mean(X)*coeff_OLS,sigma_OLS) - normcdf(lb,mean(X)*coeff_OLS,sigma_OLS);
c_initial_coeff = coeff_OLS/retain_portion;
c_initial = [c_initial_coeff; log(sigma_OLS)];

options = optimset('LargeScale','off','MaxFunEvals',10000,'Display','off');
try
    [estimator_big,log_like,exitflag,output,Gradient,Hessian] = fminunc(@(c)ML_TOBIT(c,lb,ub,X_lb,X_ub,Y_continous,X_continous,nreg,num_continous),c_initial,options);
catch
    disp('')
    disp('Oooopse, Matlab optimization toolbox is not installed, or it experienced an error.')
    disp('You may download a compatibility package on my website.')
    disp('http://www.public.iastate.edu/~hqi/toolkit')
    error('Program exits.')
end
estimator = estimator_big(1:nreg);
sigma = estimator_big(end);

if verbose
    disp(['Log likelihood ', num2str(-log_like)])
    disp('')
end

cov_Hessian = inv(Hessian);
std_c = sqrt(diag(cov_Hessian));
t_stat = estimator_big./std_c;

ME1 = (normcdf((ub-mean(X)*estimator)/sigma) - normcdf((lb-mean(X)*estimator)/sigma))*estimator;
ME2 = mean((normcdf((ub-X*estimator)/sigma) - normcdf((lb-X*estimator)/sigma)))*estimator;

result = cell(nreg+2,6);
result(1,:) = {'','Estimator','SE','z-stat','ME(avg. data)','ME(ind. avg.)'};
for m = 1:nreg
    result(m+1, 1) = {['C(',num2str(m-add_constant),')']};
    result(m+1, 2:6) = {estimator(m), std_c(m), t_stat(m), ME1(m), ME2(m)};
end
result(nreg+2, 1) = {'sigma'};
result(nreg+2, 2:6) = {sigma,std_c(m+1), t_stat(m+1), nan, nan};

disp('')
%disp(result)

end

%----------------------------------
% subfunction ML_TOBIT
%----------------------------------
function log_like = ML_TOBIT(c,lb,ub,X_lb,X_ub,Y_continous,X_continous,nreg,num_continous)

coeff = c(1:nreg);
logsigma = c(nreg+1);

if exp(logsigma) <= 0
    log_like = -10000000;
else
    resid = Y_continous -X_continous*coeff;
    RSS = resid'*resid;
    log_like_continuous = -num_continous*log(exp(logsigma)) - 0.5*RSS/exp(logsigma)^2;
    
    log_like_lower = sum(log(normcdf(lb,X_lb*coeff,exp(logsigma))));
    log_like_upper = sum(log(normcdf(-(ub-X_ub*coeff)/exp(logsigma))));
    
    log_like = log_like_continuous + log_like_lower + log_like_upper;
end
log_like = -log_like;

end