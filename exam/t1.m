%% Task 1
% Statistical Methods
% Home Exam, Spring 2014
%
% Mikael Mieskolainen

clear; close all;


%% i.)

% We have here a classical target for the one-sample Student's t-test
% We want to know if these two sample have the same mean or not. Let us
% have random variable X, which is the egg weight.
% http://fi.wikipedia.org/wiki/Studentin_t-testi

% H_0 : same mean, mu_x = mu_0
% H_1 : different mean, i.e. vitamin did work
%
% t = (x_bar - mu_0) / sqrt(var_X / n) 
n = 25;
x_bar_minus_mu_0 = 3;
var_X = 10^2;
n_dof = n - 1;  % Degrees of Freedom

% t-value
t = x_bar_minus_mu_0 / sqrt(var_X / n);

% p-value by integrating t-distribution
P_value = 1 - tcdf(t, n_dof);

fprintf('Task 1. i.) Students t-test gives a P-value: %0.3f \n', P_value);

% By choosing a significance level alpha = 0.05, the result is not
% significant enough. However, P-value is quite close, so it's better
% to carry on feeding the vitamins and collect a larger sample size.


%% ii.)

% Let's use binomial probability here. We do the "can you see the Northern
% lights" test n = 5 times. The probability to have lights is 0.8 (4/5) and
% the probability to have cloudless sky is 0.5. These physical phenomena
% are assumed to be independent, so total probability to see lights is then
% by multiplication rule: p = 0.8 x 0.5 = 0.4

P_bin = @(k, n, p) ...
      factorial(n) / (factorial(k)*factorial(n-k)) * p^k * (1-p)^(n-k);
n = 5;
p = 0.4;

% Probability to see the lights 5 nights in a row
P5 = P_bin(5, n, p);
fprintf('Task 1. ii.) Probability to see 5 nights in a row is %0.3f \n', P5);

% Probability to not to see the lights 5 nights in a row
P0 = P_bin(0, n, p);
fprintf('Task 1. ii.) Probability to not see anything is %0.3f \n', P0);


%% iii.)

% A charged particle beam consists of 95 % of pions and 5 % of kaons
% ~ Prior probabilities
P_pion = 0.95;
P_kaon = 0.05;

% P(identify as X | is a kaon)
P_kaon_kaon = 0.90;
P_pion_kaon = 0.05;
P_ndec_kaon = 0.05;

% P(identify as X | is a pion)
P_pion_pion = 0.80;
P_kaon_pion = 0.05;
P_ndec_pion = 0.15;

% Now the posteriori probabilities are:

% a.)
denom_a = (P_pion_kaon * P_kaon + P_pion_pion * P_pion);

% P(is a pion | identified as pion)
fprintf('Task 1. iii.) a.) P_pion = %0.4f \n', ...
        (P_pion_pion * P_pion) / denom_a);

% P(is a kaon | identified as pion)
fprintf('Task 1. iii.) a.) P_kaon = %0.4f \n', ...
        (P_pion_kaon * P_kaon) / denom_a);

% b.)
denom_b = (P_kaon_kaon * P_kaon + P_kaon_pion * P_pion);

% P(is a pion | identified as kaon)
fprintf('Task 1. iii.) b.) P_pion = %0.4f \n', ...
        (P_kaon_pion * P_pion) / denom_b);

% P(is a kaon | identified as kaon)
fprintf('Task 1. iii.) b.) P_kaon = %0.4f \n', ...
        (P_kaon_kaon * P_kaon) / denom_b);

% c.)
denom_c = (P_ndec_kaon * P_kaon + P_ndec_pion * P_pion);

% P(is a pion | identified as no decision)
fprintf('Task 1. iii.) c.) P_pion = %0.4f \n', ...
        (P_ndec_pion * P_pion) / denom_c);

% P(is a kaon | identified as no decision)
fprintf('Task 1. iii.) c.) P_kaon = %0.4f \n', ...
        (P_ndec_kaon * P_kaon) / denom_c);


