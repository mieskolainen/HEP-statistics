% Statistical Methods in Physics Spring 2014
% University of Helsinki
%
% Exercise 2
% Mikael Mieskolainen

close all;
clear;

fprintf('\n');

%% Task 1

% We have exponential time distribution with \lambda = 30 minutes
% (Using Poisson makes this asumption always!)
%
% The expected value of exponential distribution is 1/lambda
% 1/lambda = 1 / 2 hours = 3 sixth of an hour (10 min / 60 min = 1 of six)
% => lambda = 1/3 (the expectation value of Poisson distribution)


%% a.)

%% i.)

lambda = 1/3;
P = poisscdf(0,lambda); % analytically also trivial: exp(-lambda)
fprintf('1.) a.) i.) Probability to get away with it: %0.3f \n', P);


%% ii.)

% CDF = exp(-lambda*sum_i=0^floor(k) lambda^i / i!

% (probability of getting caught = 1 - not getting caught)
% P = 1 - 0.1 = 0.9
% =>
lambda_b = -log(0.9);

% Now scale the original ten minutes interval
t = (lambda_b / lambda) * 10;

fprintf('1.) a.) ii.) The nap could be maximally: %0.3f minutes \n \n', t);


%% b.)

%% i.)

% Poisson and Gaussian distributions

figure;

lambda = 5;      % E_poisson[X] = Var_poisson[X] = 5
mu = 5;          % Gaussian mean
sigma = sqrt(5); % Gaussian std

X = -5:0.1:15;
Y = normpdf(X,mu,sigma);
plot(X,Y); hold on;

X = 0:1:15;
Y = poisspdf(X,lambda);
stem(X,Y,'r');

xlabel('X'); ylabel('pdf(X)');
legend('Normal(5,sqrt(5))', 'Poisson(5)')


%% ii.)

% Binomial distributions
figure;

P = [0.5 0.005]; % Probability
N = [10 1000];   % Number of trials
colors = {'b','r'};

X = 0:1:15;
for i = 1:2
    Y = binopdf(X, N(i), P(i));
    stem(X,Y,colors{i}); hold on;
end

X = -5:0.1:15;
Y = normpdf(X,mu,sigma);
plot(X,Y);

legend('Bin(P = 0.5,N = 10)', ...
       'Bin(P = 0.005, N = 1000)', 'Normal(5, sqrt(5))');
xlabel('X'); ylabel('pdf(X)');


% The discrepancies are that both binomials have the same mean, but the
% second one has smaller variance. 

% 
% The second one with P = 0.005 and N = 1000 is well approximated by the
% Gaussian. 


%% iii.)

% Binomials with Poisson
figure;

colors = {'b','r+'};

X = 0:1:15;
for i = 1:2
    Y = binopdf(X, N(i), P(i));
    stem(X,Y,colors{i}); hold on;
end

stem(X, poisspdf(X,5), 'g');

legend('Bin(P = 0.5,N = 10)', 'Bin(P = 0.005, N = 1000)', 'Poisson(5)');
xlabel('X'); ylabel('pdf(X)');


%% Task 2

%% a.)

figure;

n = 10000;
r = MLCG(n);

nbins = 10;
[nout,xout] = hist(r, nbins);
bar(xout,nout,1);
xlabel('X'); ylabel('pdf(X)'); title('MLCG output');


%% i.)
% Quite uniform, not exactly due to relatively small value of n
% Lets devise some method to test uniformity. Lets calculate relative
% RMS error between expected distribution and obtained
exp = ones(nbins,1) * n / nbins;
obt = nout;
RRMS_e = sqrt(mean( ((exp(:) - obt(:)) ./ exp(:)).^2 ));
fprintf('2.) a.) i.) Relative RMS error: %0.4f \n', RRMS_e);


%% ii.)
% Mean should be: (1 - 0) / 2 = 0.5
% Variance should be: (1 - 0) / 12 = 0.0833
fprintf('2.) a.) ii.) Obtained mean (expected 0.5): %0.4f \n', mean(r));
fprintf('2.) a.) ii.) Obtained variance (expected 0.0833): %0.4f \n', var(r));


%% iii.) 


% Lets inspect the autocorrelation function, which
% should be Dirac peak if the generator does not produce
% linear correlation in the sequence

% The correlation coefficient between two preceeding values is
fprintf('2.) a.) iii.) Correlation coefficient: %0.4f \n', ...
        corr(r(2:end),r(1:end-1)));
% Which can be read from the autocorrelation sequence from also

figure;

max_delay = 50;

% NOTE, we need to substract mean here for [0,1] distributions,
% because xcorr function implementation in Matlab does not substract mean!

[c_r, lags] = xcorr(r - mean(r),max_delay,'coeff');
subplot(1,2,1);
stem(lags, c_r); axis square;
title('MLCG'); xlabel('delay'); ylabel('autocorrelation');

rr = rand(n,1);
[c_r, lags] = xcorr(rr - mean(rr),max_delay,'coeff');
subplot(1,2,2);
stem(lags, c_r); axis square;
title('rand'); xlabel('delay');

% A reasonable autocorrelation behavior for the all generators


%% b.)

figure;

RM = rand(100000,12);
Rsum = sum(RM,2);

nbins = 200;
magic = 8; % Magic normalization factor to normalize histogram
           % with respect to continuous pdf (Gaussian here)
[nout,xout] = hist(Rsum, nbins);
bar(xout, nout / sum(nout) * nbins / magic, 1);
xlabel('X'); ylabel('pdf(X)'); hold on;

title('Distribution of a sum of 12 uniformly distributed random variable');

% From CLT, the mean should be \sum_{i=1}^12 0.5 = 6
%           and variance \sum_{i=1}^12 1/12 = 1,
% lets plot a Gaussian with the same values
mu = 6; sigma = 1;
x = 0:0.1:12;
plot(x, normpdf(x,mu,sigma),'r', 'LineWidth', 1.5);


%% Task 3


%% i.)

% Description of the problem:

% Angle \theta ~ U([-\pi,pi]) (uniform)

% Lets define: \tan \theta = a/b = (X - x0) / y0
% <=> \theta = arctan((X - x0) / y0)

% We can obtain the PDF function as the derivative d\theta/dX
% f(X) = N d\theta / dX = N (d atan((X - x0)/y0)) /dX
%      = N 1 / (y0(1 + ((X - x0)/y0)^2 ))
%
% where N \in R is yet unknown normalization factor

% Now lets fix the normalization by
% \int_-\infty^\infty f(X) dX = N \pi = 1 <=> N = 1/\pi

% The cumulative distribution function is obtained by definition
% by integrating the PDF up to X

% F(X) = \int_-\infty^X f(z) dz
%      = 1/\pi \int_-\infty^X 1/(y0(1 + ((u-x0)/y0)^2) dz
%      = 1/\pi [atan((2z - 2x0) / (2y0))]_-\infty^X
%      = 1/2 + 1/\pi atan((X - x0)/y0)

% The inverse transform method works as follows. We generate random
% numbers uniformly between [0,1], and use a (inverse)
% transform function E(r), which transform these numbers to
% follow our distribution f(X).

% Lets denote uniform pdf with g(r) and its 
% cumulative distribution with G(r) = r, with r \in [0,1]

% Now we can solve the (inverse) transform function E(r) by

% F(E(r)) = G(r) = \int_-\infty^r g(z) dz = r <=>
% 1/2 + 1/\pi atan((E(r) - x0)/y0) = r <=>
% E(r) = y0(tan((r - 1/2)pi) + x0)

% Now we have all necessary to proceed.


%% ii.)

% Plot 10000 X values using the method with (x0,y0) = (0,1)

x0 = 0;
y0 = 1;

n = 10000;
r = rand(n,1);
X = y0*(tan((r - 0.5)*pi) + x0);

figure;
Xs = sort(X);
cut = 0.03; % We need to cut out outliers because this distribution is
% well-known to behave badly with the number of samples n-> infty !!
hist(Xs(cut*n:end-cut*n), 60); axis tight;
xlabel('X'); ylabel('pdf(X)'); title('Cauchy distribution');



%% iii.)

% The distribution has theoretically no expectation value
% and its variance is infinite, but it has two parameters
% x0 or (or t)  (peak position ~ mode / median)
% gamma (or s)  (full width at half maximum, FWHM)

% Okay, however we expect the central value to be 0 (the value of x0)
% and the FWHM to be 2 (twice the value of y0), i.e. gamma = 2 x y0


%% iv.)

% Fit using closed form maximum likelihood estimates for mean 
% and variance for the central part (by cutting out rest of the values)
% i.e. we use simple "robust estimation" methods here.

cut = 0.25; % Percentile to cut out of both tales,
%             to obtain FWHM boundaries
%             (http://en.wikipedia.org/wiki/Cauchy_distribution

% Truncated estimater for the central parameter, could also be obtained
% just by taking median, which is even more robust,
% i.e. the so-called fracture point is at 50%
x0_est = mean(Xs(cut*n:end-cut*n)); % Xs is sorted

% Quartile estimator for the width parameter
y0_est = Xs(end-cut*n) - Xs(cut*n);
fprintf('\n3.) iv.) x0: %0.3f (expect. 0), 2 x y0: %0.3f (expect. 2) \n', ...
        x0_est, y0_est);

% So, they agree quite ok!

% But because the exercise asked the Gaussian fit, which is simply done
% by mean and variance calculations for the truncated part. These are
% closed form maximum likelihoods estimators "fitters" for the
% Gaussian distribution, so here they are:

% Mean was calculated above (x0_est), the standard deviation is
sigma = std(Xs(cut*n:end-cut*n));

fprintf('3.) iv.) mu: %0.3f, sigma: %0.3f \n', x0_est, sigma);




