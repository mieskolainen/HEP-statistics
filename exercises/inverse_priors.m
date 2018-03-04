% Maximization of marginal likelihood
% "Process fraction fitter test"
% 
% One should extend this more, proper distribution statistics
% (multinomial, Poisson etc.), extended maximum likelihood (fitting the
% overall normalization) etc.
%
% Mikael Mieskolainen, 2014

% Class 1
mu(1) = 5;
std(1) = 1;
data1 = [std(1)*randn(1000,1) + mu(1); randn(200,1)];

% Class 2
mu(2) = 2;
std(2) = 2;
data2 = std(2)*randn(300,1) + mu(2);

% Class 3
mu(3) = 4;
std(3) = 1;
data3 = std(3)*randn(300,1) + mu(3);

X = [data1; data2; data3];

% Priors
True_priors = [length(data1) length(data2) length(data3)] / length(X)

MSE = [];
min_MSE = inf;
best_priors = [];

L = [];

s_range = 0.01:0.01:0.99;

priors_saves = zeros(length(s_range), 3);

g = 1;
for s = s_range

    Priors = [s 0.5-s/2 0.5-s/2]
    priors_saves(g,:) = Priors;
    
% Do the classification
Ps = zeros(length(X),3);
Posteriors = zeros(length(X),3);
for i = 1:length(X)

    P = zeros(1,3);
    
    P(1) = (normpdf(X(i),mu(1),std(1)) + 0.2*normpdf(X(i),0,1))/1.2*Priors(1);
    for k = 2:3
        P(k) = normpdf(X(i), mu(k), std(k)) * Priors(k);
    end
    
    % Save unnormalized
    Ps(i,:) = P;
    
    % Normalize
    Posteriors(i,:) = P / sum(P);
end

    % Calculate MSE
    MSE(g) = mean( ((mean(Posteriors) - Priors)).^2./Priors );
    if (MSE(end) < min_MSE)
       min_MSE = MSE(end);
       best_priors = Priors;
    end
    
    % Calculate log-likelihood (log for numerical reasons)
    
    % Empirical Bayes, i.e. f(x|c_j) x P(C_j), marginal likelihood, 
    % evidence maximization
    lh = 1;
    for i = 1:length(X)
        lh = lh + log(sum(Ps(i,:)));
    end
    L(g) = lh;
    
    g = g + 1;
end

True_priors
best_priors


%%

figure;
subplot(1,2,1);
plot(s_range, log(MSE)); axis square; xlabel('s'); ylabel('ln(MSE)');
subplot(1,2,2);
plot(s_range, L); axis square; xlabel('s'); ylabel('Likelihood');
%xlabel('P');
%ylabel('log_{10}(MSE)');
%axis tight;


