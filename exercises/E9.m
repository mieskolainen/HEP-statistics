% Statistical Methods in Physics Spring 2014
% University of Helsinki
%
% Exercise 9
% Mikael Mieskolainen
clear; close all;

format long;

% Poisson likelihood function for signal + background with mean = s + b
L = @(s,b,n) (s+b)^n*exp(-(s+b))/factorial(n);

% Expected number of background events (assumed to be known exactly)
b = 3.6;

% We observe the total number of events
n = 16;


%% Task 1

%% i.)

% Null hypothesis (s = 0), i.e. background only (no signal or new physics)

% We have a relation for the cumulative sum of Poisson probabilities
% Note, I changed the index to k, because n was already in use!!
% \sum_{k=0}^m P(k; v) = 1 - F_{\chi^2} (2v; n_dof)

% Okay, v = s + b
v = 0 + b;

% n_dof = 2(m+1) = 2*(n-1+1)
% (m is an auxialary variable which goes up to n - 1)
n_dof = 2*(n-1+1);

% P-value is: 1 - Poisson_CDF
P = 1 - (1 - chi2cdf(2*v, n_dof));

fprintf('i.) P-value: %0.3f x 10^(-6) \n', P*1e6);


%% ii.)

% This is a story about one-tail vs. two-tail tests

% Okay, in Particle Physics it's always understood that the 5 sigma level
% corresponds P-value of P = 1 - normcdf(5) ~= 1/3500000 ~= 0.00003 %
% i.e. we have 0.00003% on the right side of normal distribution CDF
% Check e.g. ATLAS paper:
% http://cdsweb.cern.ch/record/1421964/files/science.pdf?version=1
% where they note 1.4 % = 2.2 sigmas,
% i.e. sigma = norminv(1 - 0.014) = 2.2

% The Poisson hypothesis test above was _one-tailed_ test, because "if the
% test statistic is always positive (or zero), only the one-tailed test is
% generally applicable", i.e. we have counts in Poisson measurement
% (which are always >= 0)
% http://en.wikipedia.org/wiki/One-_and_two-tailed_tests

% In general, chi^2 tests are one-tailed tests, because chi^2 distribution
% is one tailed. Also, it makes physical sense here because we assumed that
% we know background exactly, so the devition should be upwards only
% - the possible signal! Thus, so by using inverse of the Normal CDF
% (norminv) we get:
% >> THIS IS CORRECT HERE
one_tailed_sigmas = norminv(1-P);

% The two-tailed tests of 5 sigmas corresponds to change
% P = 1 - (normcdf(5) - normcdf(-5)) ~= 1/1750000 ~= 0.00006 % change
% i.e. we have 0.00003% on the left side and 0.00003% on the right side of
% the normal distribution CDF. Using this convention, we get sigmas by
% dividing our one-tailed P-value by 2 (this division is valid because
% Normal distribution is symmetric distribution, and the division
% corresponds to substraction of the left tail mass in the Normal CDF)
% >> THIS IS "WRONG" HERE
two_tailed_sigmas = norminv(1-P/2);

fprintf('ii.) One-tailed (the correct): %0.3f x sigma, two-tailed: %0.3f x sigma \n', ...
        one_tailed_sigmas, two_tailed_sigmas);

% Summa summarum,
% the two tailed meaning of 5 sigma limit for discovery is weaker!

% Anyway, no discovery here (i.e. we don't reject the null hypothesis of
% background only, because we have this very strict criteria of 5 sigmas).


%% iii.)

% The notation in the exercise sheet is still non-correct,
% (at least the left side of the ratio formula), I correct it below:

% Now assume we have seen a signal of a new physics phenomena
% Calculate now the confidence interval for the number of signal
% events at 68.3 (1 sigma) confidence level.
 
% We calculate CL interval by variating s in the likelihood ratio:
%
% \lambda(s; b, n) = L(s, b; n) / L(s_hat, b; n)

% b is exact here, i.e. that doesn't change and n is the number of
% measurements, the ratio is a function of s only (n is not varied).

% s_hat is estimated by maximizing likelihood of the data under known
% background i.e. d ln L/ds = 0 => which gives (Mathematica) 
% -exp(-b-s)(b+s)^(n-1)(b-n+s) = 0 =>
s_hat = n - b;
% Very intuitive results for the ML estimate of number of signal events.

% Note! Error in exercise sheet which states that we should differentiate
% with respect to b. No, b was exact, we need to differentiate with
% respect to s to get s_hat. Well, luckily it gives the same answer
% just because of symmetry in the formula between s and b...

% Likelihood ratio function (simplified by using s_hat = n - b) =>
% ( (s+b)^n*exp(-(s+b))/n! ) / ( (s_hat+b)^n*exp(-(s_hat+b))/n! ) =>
lambda = @(s, b, n) (s+b)^n / (n^n) * exp(n-s-b);

% Now variate s ( b and n are fixed as my notation lambda(s; b, n) states)
s_range = 6:0.001:18;
minus2lnLambda = zeros(length(s_range),1);
for i = 1:length(s_range)
    minus2lnLambda(i) = -2*log(lambda(s_range(i), b, n));
end

% Find minimum position (=corresponds to s_hat)
[min_value, min_ind] = min(minus2lnLambda);

% Find CL-limits where as points >= -2lnLambda + 1
% (corresponds to the point chi2cdf(1, ndof) with ndof = 1)
% (basically the usual lnL - 0.5)
% For more info, see:
% http://arxiv.org/pdf/physics/0403059v5.pdf
CL_low = 0;
for i = min_ind:-1:1
    if (minus2lnLambda(i) >= min_value + 1)
        CL_low = s_range(i);
        break;
    end
end
CL_high = 0;
for i = min_ind:length(s_range)
    if (minus2lnLambda(i) >= min_value + 1)
        CL_high = s_range(i);
        break;
    end
end

% Plot it out
figure;
plot(s_range, minus2lnLambda); hold on;
plot(s_hat*ones(2,1), ...
     [min(minus2lnLambda) max(minus2lnLambda)], 'k-');
plot(CL_low*ones(2,1), ...
     [min(minus2lnLambda) max(minus2lnLambda)], 'r--');
plot(CL_high*ones(2,1), ...
     [min(minus2lnLambda) max(minus2lnLambda)], 'r--'); 
plot(s_range, ones(size(s_range)), 'k--');
xlabel('Signal rate s');
ylabel('Profile likelihood -2ln \Lambda(s; b,n)');
title('Likelihood ratio test \Lambda, ML estimate of s and its 1\sigma CL-bounds (red)');
axis tight;

fprintf('iii.) s_hat = %0.2f with 1 sigma CL bounds: s = [%0.2f, %0.2f] \n',...
        s_hat, CL_low, CL_high);


%% iv.)

N = 1000; % Number of pseudoexperiments

% Random background event numbers with mean b and std below
std_b = 1.2;
b_r = randn(N,1)*std_b + b;

% Now variate s
s_range = 6:0.01:18;
minus2lnLambda = zeros(N, length(s_range));

% Generate "pseudoexperiments"
for i = 1:N
    for j = 1:length(s_range)
        minus2lnLambda(i,j) = -2*log(lambda(s_range(j), b_r(i), n));
    end
end

% Find minimum at each s value
[min_values, min_inds] = min(minus2lnLambda, [], 1);

% Now find the mean of minimum likelihood ratios
min2lnLambda_min = mean(min_values);


%% 

% Calculate 2D histogram of the pseudoexperiements for each s
bins = min(minus2lnLambda(:)):0.1:quantile(minus2lnLambda(:), 0.999);
histM = zeros(length(bins), length(s_range));

for j = 1:length(s_range)
    [histM(:,j), ~] = hist(minus2lnLambda(:,j), bins);
end
histM = histM / sum(histM(:)); % Normalize 2D

% Plot lnL surface
figure;
surf(s_range, bins, histM, 'EdgeColor', 'none'); colorbar;
axis tight;
xlabel('Signal rate s');
ylabel('-2ln \Lambda');
zlabel('pdf(s, 2-ln \Lambda)');
view([65,40]);
title('Part iv.)');


%%

% Calculate the number of pseudoexperiments for each s which cover 68.3%
frac = zeros(length(s_range),1);
for j = 1:length(s_range)
    for i = 1:N
        if (minus2lnLambda(i,j) <= min2lnLambda_min + 1)
            frac(j) = frac(j) + 1;
        end
    end
end
frac = frac / N;

% Find the CL-bounds (0.683)
CL_low  = s_range(find(frac >= 0.683, 1, 'first'));
CL_high = s_range(find(frac >= 0.683, 1, 'last'));

figure;
plot(s_range, frac); hold on;
plot(s_range, ones(size(s_range))*0.683, 'k--');

% Plot it out
plot(s_hat*ones(2,1), [0 1], 'k-');
plot(CL_low*ones(2,1), [0 1], 'r--');
plot(CL_high*ones(2,1), [0 1], 'r--'); 
xlabel('Signal rate s');
ylabel('Fraction of passed pseudoexperiments');
axis tight;


fprintf('iv.) New 1 sigma CL bounds: s = [%0.2f, %0.2f] \n',...
        CL_low, CL_high);

% Slightly smaller CL-region than in iii.)

