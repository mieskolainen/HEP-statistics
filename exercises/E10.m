% Statistical Methods in Physics Spring 2014
% University of Helsinki
%
% Exercise 10
% Mikael Mieskolainen

% Histogram unfolding example using simple Tikhonov regularization
% approach. Method does not enforce explicitly non-negativity of
% the estimated histogram, but works here for demonstration purposes.

clear; close all;

% ------------------------------------------------------------------------
% Signal model:   m = a + s, where s is noise, m is measurement, a is
%                            original signal
% 
% Folded mapping: M = W*A, where M and A are histogram vectors and W is
%                          the system (smearing) matrix
% 
% Technical note! Signal model above is not per se "smearing", but just
% additive noise model. However, it can be interpreted as
% smearing ~ convolution of the corresponding probability distributions:
% 
% "The probability distribution of the sum of two or more
% independent random variables is the convolution of their individual
% distributions."
% - en.wikipedia.org/wiki/Convolution_of_probability_distributions
%
% Thus:     pdf(m) = pdf(a + s)
%       <=> pdf(m) = pdf(a) (*) pdf(s), where (*) denotes convolution
%
% And this makes this problem a deconvolution problem, because in this
% case we want to solve the pdf(a), not the actual signal vector a.
% That would be a denoising problem!
% ------------------------------------------------------------------------


%% i.) The "Monte Carlo simulation" part

% ------------------------------------------------------------------------
% 1. Simulate measurements

a = randn(1000000,1)*2.5 - 0.2;   % Create original "true" signal a
s = randn(1000000,1)*0.4 + 2.0;   % Create additive noise signal s
m = a + s;                        % Create simulated measurement m

% ------------------------------------------------------------------------
% 2. Calculate histograms (distributions), A and M

bins = -10:1/3:10;

% Now create histograms, could be normalized if wanted by: X / sum(X(:))
A = hist(a, bins)';
M = hist(m, bins)';

% ------------------------------------------------------------------------
% 3. Calculate system characterization (smearing) matrix W

% Two-dimensional histogram
W = hist3([m a], {bins, bins});

% Normalize each row
for i = 1:size(W,1)
    W(i,:) = W(i,:) / sum(W(i,:));
end;

figure;
imagesc(W);
xlabel('j');
ylabel('i');
title('Smearing matrix W: A -> M, where [A]_{ij} = P(i | j)');
colorbar;
axis square;

figure;
stem(bins, A, 'b'); hold on;
stem(bins, M, 'r');

legend('True, A', 'Smeared, M'); legend('boxoff');
xlabel('Signal value');
ylabel('Histogram count (#)');
title('Part i.)');


%% ii.) The "measurement" part

a = randn(100000,1)*2.0 + 0.0;
s = randn(100000,1)*0.4 + 2.0;
m = a + s; % smearing
A = hist(a, bins)';
M = hist(m, bins)';

figure;
stem(bins, A, 'b'); hold on;
stem(bins, M, 'r');

legend('True, A', 'Smeared, M'); legend('boxoff');
xlabel('Signal value');
ylabel('Histogram count (#)');
title('Part ii.)');


%% iii.) The "unfolding the measurement" part

% Solve the unfolding problem inv(W) : M -> A

% Solve "unfolding" by classic Tikhonov regularization
% M = W*A  <=>  A_hat = pinv(W)*M,
% where pinv() is the pseudoinverse operator with regularization term

lambda = 0.02; % Regularization strength, play with this, might fail

% Construct second order derivative (finite difference) matrix,
% i.e. curvature, to induce explicit smoothness into the solution
L = zeros(length(M));
for i = 1:size(L,1)-2
   L(i,i:i+2) = [1 -2 1];
end
%L = eye(length(M));

% Signal covariance matrix
V = eye(length(M));

% http://www.desy.de/~blobel/school_march10.pdf
% \ operator does not work here, use standard inv()
A_hat = inv(W'*V*W + lambda*(L'*L)) * W'*V * M; % <- MAGIC HAPPENS!

% RMSE between distributions
fprintf('iii.) RMSE before unfolding: %0.1f, after unfolding: %0.1f \n', ...
        sqrt(mean((M - A).^2)), sqrt(mean((A_hat - A).^2)));

% Plot all
figure;
stem(bins, A, 'b'); hold on;
stem(bins, M, 'r');
stem(bins, A_hat, 'g.-');
axis tight;

legend('True, A', 'Smeared, M', 'Unfolded, A_{hat}'); legend('boxoff');
xlabel('Signal value'); ylabel('Histogram count (#)'); title('Part iii.)');


% Comments: I avoided using the same mapping W in the inverse simulation,
% because then we would have done the so-called "inverse crime" simulation,
% i.e. that's why I used different noise realization in ii.) than in
% i.) where we estimated the mapping W. Makes the problem more realistic
% and the inversion more unstable.

