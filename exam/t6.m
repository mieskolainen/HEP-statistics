%% Task 6
% Statistical Methods
% Home Exam, Spring 2014
%
% Mikael Mieskolainen

clear; close all;

warning off;

% Time points (ns)
global t; % Global because of the optimization routine
t = [300 350 400 450 500 550 600 650 700 750]';

% Number (N) of decays in 10 ns windows at the given time points in t
N = [803000 581083 429666 320016 242783 188487 150737 124103 105397 92748]';

% Decay rates, R = dN/dT ~ deltaN /deltaT when deltaT -> 0
global R; % Global because of the optimization routine
R = N / 10; % 10 ns window,
% i.e. we approximate dN/dT using a finite difference deltaN / deltaT,
% i.e. secant approximation


%% i.)

% Let's do a non-linear least squares fit, because all the change of
% variable tricks and (log) transforms together with linear fit 
% fail in this case
% (if we want to keep lambda, A and B as a seperate variables)

% Parameter vector (lambda, A, B) initial estimate,
% this needs to be reasonably good
theta0 = [-0.006 5.5e5 5000];

% Cost function
costf = @(theta) R - ( theta(2)*exp(-theta(1)*t) + theta(3) ); 

% Non-linear Least Squares, basically it does vector Taylor's polynomial
% expansion => Jacobian matrix linearization of the usual normal
% equations of the linear least squares, which are iterated into the
% direction of the gradient of the least squares cost.
% Guaranteed to converge to a local minimum only in this case,
% because the problem is non-convex, thus the initial estimate must be good.

% Jacobian matrix is approximated by the routine with finite differences,
% which is good in this case because our function is smooth.
[theta_hat,resnorm,residual,~,~,~,jacobian] = lsqnonlin(costf, theta0);
jacobian = full(jacobian); % Make it non-sparse datatype

% Parameter estimates
lambda = theta_hat(1);
A = theta_hat(2);
B = theta_hat(3);

% One of the reason to use non-linear least squares here with Jacobian
% matrix information is the easy uncertainty estimation based on
% RCF-theory.

% Uncertainty estimates by linearization (G. Gowan), i.e. using the
% Jacobian matrix (Fisher's) information evaluated at the parameter
% estimate point (lambda_hat, A_hat, B_hat) is

% Parameter covariance matrix estimate
% Note that we use residual as error estimates here, as often done.
C = mean(residual.^2)*inv(jacobian'*jacobian)';
std_lambda = sqrt(C(1,1));
std_A = sqrt(C(2,2));
std_B = sqrt(C(3,3));

% Parameter estimates and their uncertainties
fprintf('Task 6. i.) lambda = %0.7f +- %0.7f \n', lambda, std_lambda);
fprintf('Task 6. i.) A = %0.0f +- %0.4f \n', A, std_A);
fprintf('Task 6. i.) B = %0.2f +- %0.5f \n', B, std_B);

% However, these parameter estimates are based on the local optimum
% estimate of the parameter vector, and their uncertainties on 
% the LS-cost function's curvature around it. Actually we should use first
% some global optimization technique like
% simulated annealing to find the best initial estimate for the non-linear
% least squares iteration, i.e. the estimated uncertainties do not take
% per se into account the cost functions multiple minimas (that's why they
% are objectively speaking relatively small).

% Plot data
figure;
plot(t, R, 'r.', 'MarkerSize', 14); hold on;

% Plot fit
t_range = 250:800;
plot(t_range, A*exp(-lambda*t_range) + B, 'k-');
axis tight;

xlabel('t (ns)'); ylabel('R(t)');
legend('Data', 'Non-Linear LS-Fit')


%% ii.)

% Lifetime, as given in the exam paper as a function of decay constant,
% NOTE! I think this is actually half-life, not mean lifetime which is
% by definition tau = 1 / lambda. However, let's calculate both.
tau = log(2) / lambda;
tau_alternative = 1 / lambda;

% Wikipedia http://en.wikipedia.org/wiki/Exponential_decay

% The uncertainty on tau can be obtained easily by error propagation:
% d/d lambda ( ln(2) / lambda ) = -ln(2) / lambda^2

% Error propagation formula gives the uncertainty on tau as
sigma_tau = sqrt((-log(2) / (lambda^2))^2 * std_lambda^2);
sigma_tau_alternative = sqrt((- 1 / (lambda^2))^2 * std_lambda^2);

fprintf('Task 6. ii) Lifetime (= ln 2/lambda) is tau = %0.1f +- %0.3f (ns) \n', ...
        tau, sigma_tau);
fprintf('Task 6. ii) Lifetime (= 1/lambda) is tau = %0.1f +- %0.3f (ns) \n', ...
        tau_alternative, sigma_tau_alternative);
    
% The 1/lambda way to calculate lifetime seems to give a value which is
% much closer to the Wikipedia mean lifetime (138 ns) of the orthopositronium

