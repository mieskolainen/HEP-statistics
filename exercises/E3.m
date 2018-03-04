% Statistical Methods in Physics Spring 2014
% University of Helsinki
%
% Exercise 3
% Mikael Mieskolainen
clear; close all;


%% TASK 1

% Photomultiplier dynode cascade
%           n1    n2           n6
% e^- -> () -> () -> ... -> () ->  

% Poisson variable n_i ~ Poiss(v_i)

N = 6;      % Number of dynodes
v = 3.2;    % Expectation value for each i-th dynode
S = 10000;   % Number of simulations with one incoming photoelectron (PE)

% Simulation results matrix
R = zeros(S, N);

%% i.)

% Loop over simulations
for k = 1:S
    
    n = zeros(N,1);     % Dynode outputs (# electrons)
    n(1) = poissrnd(v); % The output of the first dynode
    
    % Loop over the rest of the dynodes i = 2,...,N
    % The no output case, i.e. n(d-1) = 0 is handled automatically
    for d = 2:N
        n(d) = sum(poissrnd(v, n(d-1), 1)); % Sum over each amplifications
    end
    
    % Save the resulting amplification cascade
    R(k,:) = n;
end

%% ii.)

% The number of electrons expectation value is theoretically
v_out = prod(ones(N,1)*v);

% Simulation value
v_out_hat  = mean(R(:,N)); 
sigma2_hat = var(R(:,N));

fprintf('T1 ii.) Expected mean: %0.2f, Obtained Mean: %0.2f, Var: %0.2f \n', ...
        v_out, v_out_hat, sigma2_hat);

% Histogram of n_out (n_N)
hist(R(:,N), 100);
xlabel('n_{out}');
ylabel('Count (#)');
axis tight;
title('PMT output');

%% iii.)

% The fraction of simulations with n_out = 0
f = sum(R(:,N) == 0) / S;

fprintf('T1 iii.) Fraction with n_out = 0 is: f = %0.3f \n', f);


%% iv.)

P = poisspdf(0:10, v);

% A. First gives 0
% B. First gives 1 second gives 0
% C. First gives 1 second gives 1 third gives 0
P_hat = P(1) + P(2)*P(1) + P(2)*P(2)*P(1);

fprintf('T1 iv.) Analytic estimate with 3 largest contributions: %0.3f \n', ...
       P_hat);


%% TASK 2

% Muons
mu_mu = 0;
sigma_mu = 1;

% Electrons
mu_e = 3;
sigma_e = 1;

% Lets plot likelihood functions (pdfs)
figure;
X = -6:0.1:8;
plot(X, normpdf(X, mu_mu, sigma_mu), 'r'); hold on;
plot(X, normpdf(X, mu_e, sigma_e), 'b');
xlabel('t'); ylabel('f(t | class)'); title('Class likelihoods (pdfs)');

% Draw selection criteria t < 1 for muons
t = 1;
line([t t], [0 0.5], 'Color',[0 0 0]);
legend('Muons', 'Electrons', 't < 1');

% The selection efficiency with t < 1
P = normcdf([-inf t], mu_mu, sigma_mu);
fprintf('T2 i.) Muon selection efficiency: %0.2f \n', P(2));


%% ii.)

% Probability that electron will be accepted as muon with t < 1
t = 1;
P = normcdf([-inf t], mu_e, sigma_e);
fprintf('T2 ii.) Probability that electron accept. as a muon: %0.2f \n', ...
        P(2));


%% iii.)

% Priors
pi_mu = 0.05;
pi_e = 0.95;

% Purity is defined as: N_signal / (N_signal + N_background)

% int_-inf^t prior_mu*pdf(t | muon) / ...
%        int_-inf^t prior_e*pdf(t | electron) + prior_mu * pdf(t | muon)

t = 1;
P = (pi_mu * normcdf(t, mu_mu, sigma_mu)) / ...
    (pi_e * normcdf(t, mu_e, sigma_e) + pi_mu*normcdf(t, mu_mu, sigma_mu));

fprintf('T2 iii.) Purity of the muon sample with t < 1 is: %0.2f \n', P);


%% iv.)

% Wanted purity for muons
purity = 0.8;

t = -1:1e-4:1;
P = zeros(size(t));

% Loop over possible t-values
for i = 1:length(t)
    P(i) = (pi_mu * normcdf(t(i), mu_mu, sigma_mu)) / ...
(pi_e*normcdf(t(i), mu_e, sigma_e) + pi_mu*normcdf(t(i), mu_mu, sigma_mu));
end

% Find the corresponding index
[~,ind] = min(abs(P - purity));
t_cut = t(ind);

fprintf('T2 iv.) Muon purity 0.8 obtained with: t < %0.4f \n', t_cut);

figure;

% Posterior distributions
P_mu = (normpdf(X, mu_mu, sigma_mu)*pi_mu) ./ ...
       (normpdf(X, mu_mu, sigma_mu)*pi_mu + normpdf(X, mu_e, sigma_e)*pi_e) ;

P_e = (normpdf(X, mu_e, sigma_e)*pi_e) ./ ...
       (normpdf(X, mu_mu, sigma_mu)*pi_mu + normpdf(X, mu_e, sigma_e)*pi_e) ;

plot(X, P_mu, 'r'); hold on;
plot(X, P_e, 'b');
xlabel('t'); ylabel('Posterior (class | t)');
title('Class posteriors and selection criteria');

% Draw selection criteria t < 1 for muons
t = 1;
line([t t], [0 1], 'linestyle', '-', 'Color',[0 0 0]);

% Draw Bayes optimal selection criteria
line([0.5 0.5], [0 1], 'linestyle', '--', 'Color',[0.25 0.25 0.25]);

% Draw selection criteria as wanted in iv.)
line([t_cut t_cut], [0 1], 'Color', [0.75 0.75 0.75]);

legend('Muons, \pi_{\mu} = 0.05', 'Electrons, \pi_e = 0.95', ...
       't < 1', 'Bayes MAP criteria', 'Muon purity 0.80');

% Show the Mathematica solution
figure;
I = imread('erfc.png');
imshow(imrotate(I,90), []);


%% EXTRA

% Lets plot purity and efficiency for muon and electron as a function of t_cut

t = -6:0.05:6;
P_mu = zeros(size(t));
P_e = zeros(size(t));

% Loop over possible t-values
for i = 1:length(t)
    P_mu(i) = (pi_mu * normcdf(t(i), mu_mu, sigma_mu)) / ...
(pi_e*normcdf(t(i), mu_e, sigma_e) + pi_mu*normcdf(t(i), mu_mu, sigma_mu));

    P_e(i) = (pi_e * normcdf(t(i), mu_e, sigma_e)) / ...
(pi_e*normcdf(t(i), mu_e, sigma_e) + pi_mu*normcdf(t(i), mu_mu, sigma_mu));
end

figure;

% Purities
plot(t, P_mu, 'r-'); hold on;
plot(t, P_e, 'b-');

% Efficiencies, for the electron (1 - CDF, because cut is to other
% direction)
eff_mu = normcdf(t, mu_mu, sigma_mu); 
eff_e = 1 - normcdf(t, mu_e, sigma_e);

plot(t, eff_mu, 'r--');
plot(t, eff_e, 'b--');

% Purity times efficiency
plot(t, P_mu.*eff_mu, 'r-.');
plot(t, P_e.*eff_e, 'b-.');

% Bayes rate (Bayes error = 1 - this)
plot(t, eff_mu*pi_mu + eff_e*pi_e, 'k--');

xlabel('t_{cut}');
ylabel('value');
legend('Muon purity', 'Electron purity', 'Muon efficiency', ...
       'Electron efficiency', 'Muon pur x eff', 'Electron pur x eff', ...
       'Bayes rate', 'Location', 'SouthWest');
   legend('boxoff');
title('Selection (cut) rule: t < t_{cut}, if true = muon, otherwise electron');






