% Statistical Methods in Physics Spring 2014
% University of Helsinki
%
% Exercise 5
% Mikael Mieskolainen

clear;
close all;


%% Task 1

% Read in data for theory 1
fid = fopen('./data/theory_1.dat');
T1 = fscanf(fid, '%f %f %f', [3 inf])';
t1 = T1(2:end, 1);  % Skip the first row
fclose(fid);

% Read in data for theory 2
fid = fopen('./data/theory_2.dat');
T2 = fscanf(fid, '%f %f %f', [3 inf])';
t2 = T2(2:end, 1);  % Skip the first row
fclose(fid);

% Read in data (actually theory 1 with noise, original data file missing...)
fid = fopen('./data/theory_1.dat');
D = fscanf(fid, '%f %f %f', [3 inf])';
d = D(2:end, 1);    % Skip the first row

% Add some Gaussian with truncation (just for the show)
d = d + randn(size(d))*1.0; d(d < 0) = 0;

fclose(fid);


N = length(d); % Number of Degrees of Freedom, here it's the number of bins

% Chi^2, i.e. test theory vs. data
chi2_1 = sum((d - t1).^2 ./ t1);
chi2_2 = sum((d - t2).^2 ./ t2);

% P-value is int_x^2^inf chi^2-distribution, thus 1 - cdf
p_1 = 1 - chi2cdf(chi2_1, N);
p_2 = 1 - chi2cdf(chi2_2, N);

fprintf('Task 1.) chi^2, P-value: theory_1: %0.3f, %0.3f, theory_2: %0.3f, %0.3f \n', ...
        chi2_1, p_1, chi2_2, p_2);

% Based on P-values, we reject theory_2 with p-value much less than 0.05,
% theory 1 is compatible with data with high p-value.


%% Task 2

% Calculate pseudo-experiment distribution for the theory

k = 1e4;        % The number of pseudoexperiments
N_bins = 300;   % The number of bins in histograms


%% Theory 1

chi2 = zeros(k,1);

% Generate a pseudoexperiment and calculate chi^2 for each
for i = 1:k
    ps = poissrnd(t1);                  
    chi2(i) = sum((ps - t1).^2 ./ t1);
end

[n,x1] = hist(chi2, N_bins);
pdf1 = n / k;
cdf1 = cumsum(pdf1);

figure;
set(gcf,'visible','off');

subplot(2,2,1);
plot(x1, pdf1); hold on;
plot([chi2_1 chi2_1], [0 0.05], '--k'); axis tight;
xlabel('\chi^2'); ylabel('pdf(\chi^2)');
legend('Pseudo', 'Data \chi^2'); legend('boxoff');
title('Theory 1');

% Cumulative plot
subplot(2,2,3);
plot(x1, cdf1); hold on;
plot([chi2_1 chi2_1], [0 1], '--k'); axis tight;
xlabel('\chi^2'); ylabel('cdf(\chi^2)');
legend('Pseudo', 'Data \chi^2'); legend('boxoff');


%% Theory 2

chi2 = zeros(k,1);

% Generate a pseudoexperiment and calculate chi^2 for each
for i = 1:k
    ps = poissrnd(t2);                  
    chi2(i) = sum((ps - t2).^2 ./ t2);
end

[n2,x2] = hist(chi2, N_bins);
pdf2 = n / k;
cdf2 = cumsum(pdf2);

set(gcf,'visible','on');
subplot(2,2,2);
plot(x2, pdf2); hold on;
plot([chi2_2 chi2_2], [0 0.05], '--k'); axis tight;
xlabel('\chi^2'); ylabel('pdf(\chi^2)');
legend('Pseudo', 'Data \chi^2'); legend('boxoff');
title('Theory 2');

% Cumulative plot
subplot(2,2,4);
plot(x2, cdf2); hold on;
plot([chi2_2 chi2_2], [0 1], '--k'); axis tight;
xlabel('\chi^2'); ylabel('cdf(\chi^2)');
legend('Pseudo', 'Data \chi^2'); legend('boxoff');


%% P-values

% Find the p-values by calculating the fraction of pseudoexperiments
% giving larger chi^2 than the one in the direct comparison
% between data vs. theory
[~, ind1] = min(abs(chi2_1 - x1));
[~, ind2] = min(abs(chi2_2 - x2));
p1 = 1 - cdf1(ind1);
p2 = 1 - cdf2(ind2);

fprintf('Task 2.) P-value: theory_1: %0.3f, theory_2: %0.3f \n', ...
        p1, p2);

% Based on P-values, we reject theory_2 with p-value less than alpha = 0.05,
% theory 1 is compatible with data with high p-value. Same conclusions as
% in Task 1.

