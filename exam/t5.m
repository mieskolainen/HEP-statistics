%% Task 5
% Statistical Methods
% Home Exam, Spring 2014
%
% Mikael Mieskolainen

clear; close all;

% Read in data
fid = fopen('./data/real_mass.dat');
D = fscanf(fid, '%f %f %f', [3 inf])';
d = D(:, 3);        % Counts in the third column
fclose(fid);

% Read in data for theory 1
fid = fopen('./data/MC1_mass.dat');
T1 = fscanf(fid, '%f %f %f', [3 inf])';
t1 = T1(:, 3);      % Counts in the third column
fclose(fid);

% Read in data for theory 2
fid = fopen('./data/MC2_mass.dat');
T2 = fscanf(fid, '%f %f %f', [3 inf])';
t2 = T2(:, 3);  % Counts in the third column
fclose(fid);


%% i.)

% Number of Degrees of Freedom, here it's the number of bins
N = length(d); 

% Pearson's Chi^2 test, i.e. test theory vs. data
chi2_1 = sum((d - t1).^2 ./ (t1 + eps)); % + eps so we don't divide by 0
chi2_2 = sum((d - t2).^2 ./ (t2 + eps));

% P-value = int_x^2^inf pdf(x | chi^2), thus 1 - cdf(x | chi^2)
p_1 = 1 - chi2cdf(chi2_1, N);
p_2 = 1 - chi2cdf(chi2_2, N);

fprintf('Task 5. i.) MC1: (chi^2, P-value) = %0.1f, %0.3f, \n MC2: (chi^2, P-value) = %0.2f, %0.3f \n', ...
        chi2_1, p_1, chi2_2, p_2);

% Based on P-values, we reject MC1 only with P-value much
% less than 0.05. However, MC2 only scenario is somewhat
% compatible with the data, given that the P-value is slightly
% larger than 0.05.
% 
% However, one must remember that the significance
% level alpha = 0.05 as a cut off is somewhat arbitrary, just a convention.
% Summa summarum, MC2 fits data much better than MC1, but not extremely
% well.


%% ii.)


% Calculate pseudo-experiments

k = 20000;      % The number of pseudoexperiments
N_bins = 300;   % The number of bins in histograms


% MC1 pseudoexperiments -->

%{
chi2ps1 = zeros(k,1);
% Generate a pseudoexperiment and calculate chi^2 for each
for i = 1:k
    
    % Pseudo-experiment is now our "data"
    ps = poissrnd(t1);                  
    chi2ps1(i) = sum((ps - t1).^2 ./ (t1 + eps));
end

save('./matlabdata/T5_2_1.mat', 'chi2ps1');
%}
load('./matlabdata/T5_2_1.mat');

[n1,x1] = hist(chi2ps1, N_bins);
pdf1 = n1 / sum(n1);
cdf1 = cumsum(pdf1);

figure;
set(gcf,'visible','off');

subplot(2,2,1);
semilogx(x1, pdf1); hold on;
semilogx([chi2_1 chi2_1], [0 max(pdf1)], '--k');
xlabel('\chi^2'); ylabel('pdf(\chi^2)');
legend('Pseudo', 'Data \chi^2'); legend('boxoff');
title('MC 1');

% Cumulative plot
subplot(2,2,3);
semilogx(x1, cdf1); hold on;
semilogx([chi2_1 chi2_1], [0 1], '--k');
xlabel('\chi^2'); ylabel('cdf(\chi^2)');
legend('Pseudo', 'Data \chi^2'); legend('boxoff');


% MC2 pseudoexperiments -->

%{
chi2ps2 = zeros(k,1);
% Generate a pseudoexperiment and calculate chi^2 for each
for i = 1:k
    
    % Pseudo-experiment is now our "data"
    ps = poissrnd(t2);                  
    chi2ps2(i) = sum((ps - t2).^2 ./ (t2 + eps));
end
save('./matlabdata/T5_2_2.mat', 'chi2ps2');
%}
load('./matlabdata/T5_2_2.mat');

[n2,x2] = hist(chi2ps2, N_bins);
pdf2 = n2 / sum(n2);
cdf2 = cumsum(pdf2);

set(gcf,'visible','on');
subplot(2,2,2);
plot(x2, pdf2); hold on;
plot([chi2_2 chi2_2], [0 max(pdf2)], '--k'); axis tight;
xlabel('\chi^2'); ylabel('pdf(\chi^2)');
legend('Pseudo', 'Data \chi^2'); legend('boxoff');
title('MC 2');

% Cumulative plot
subplot(2,2,4);
plot(x2, cdf2); hold on;
plot([chi2_2 chi2_2], [0 1], '--k'); axis tight;
xlabel('\chi^2'); ylabel('cdf(\chi^2)');
legend('Pseudo', 'Data \chi^2'); legend('boxoff');


% Pseudoexperiment based P-values -->

% Find the p-values by calculating the fraction of pseudoexperiments
% giving larger chi^2 than the one in the direct comparison
% between data vs. theory
[~, ind1] = min(abs(chi2_1 - x1));
[~, ind2] = min(abs(chi2_2 - x2));
p1 = 1 - cdf1(ind1);
p2 = 1 - cdf2(ind2);

fprintf('Task 5. ii.) P-value: MC1: %0.3f, MC2: %0.3f \n', ...
        p1, p2);

% Based on P-values, we reject again MC1 with p-value much less than
% alpha = 0.05. MC2 has a P-value which is barely over 0.05.
% So quantititative results are almost the same as in i.). MC2 is the best
% model for our data, but not with a very good fit in chi2 sense.

% Results are not identical however, because by using the theoretical 
% chi^2 CDF (in i.) in a case of small event sample sizes might not
% give realistic results. The whole idea of pseudoexperiments is to see
% how (Poissonian) uncertainty in the bin counts affects our
% P-value estimates.


%% iii.)

% Find the optimal alpha by minimizing chi^2

% The range of constant a (wound by a priori / posteriori inspection)
a_range = linspace(0.05, 0.15, 100);

% Do a brute force scan
chi2_a = zeros(size(a_range));
for i = 1:length(a_range)
    a = a_range(i);
    chi2_a(i) = sum( (d - (a*t1 + (1-a)*t2)).^2 ./ (a*t1 + (1-a)*t2) );
end

% Find optimal value for a
[chi2a_min, min_ind] = min(chi2_a);
a_best = a_range(min_ind);

% P-value
P_value_direct_a = 1 - chi2cdf(chi2a_min, N);

fprintf('Task 5. iii.) Direct: a = %0.3f with (chi^2, P-value) = (%0.2f, %0.3f)\n', ...
        a_best, chi2a_min, P_value_direct_a);

% Do the same as above now with pseudoexperiments

%{
k = 50000;

chi2_a_ps = zeros(k, length(a_range));
for i = 1:k
    % Two independent backgrounds
    ps1 = poissrnd(t1);
    ps2 = poissrnd(t2);    
    
    for j = 1:length(a_range)
        a = a_range(j);
        
        % Now this pseudoexperiment is our "data"
        ps = a*ps1 + (1-a)*ps2;
        chi2_a_ps(i,j) = ...
            sum( (ps - (a*t1 + (1-a)*t2)).^2 ./ ((a*t1 + (1-a)*t2)) );
    end
end

save('./matlabdata/T5_3.mat', 'chi2_a_ps');
%}
load('./matlabdata/T5_3.mat');

P_values = zeros(length(a_range), 1);
for j = 1:length(a_range)
    
    % Now calculate pseudoexperiment based chi^2 cdf
    [n_j,x_j] = hist(chi2_a_ps(:,j), N_bins);
    pdf_j = n_j / sum(n_j);
    cdf_j = cumsum(pdf_j);
    
    % Now using the chi2 value from the _direct part_ and the generated
    % chi^2 cdf distribution above to calculate more honest P-value
    [~, ind] = min( abs(x_j - chi2_a(j) ) ); % Find the point at CDF
    P_values(j) = 1 - cdf_j(ind);
end

% Find maximum of pseudoexperiment based P-values,
% and this maximum corresponds to an optimal value of a
% under pseudoexperiments
[~, opt_ind] = max(P_values);
fprintf('Task 5. iii.) PSE: a = %0.3f with (chi^2, P-value) = (%0.2f, %0.3f) \n', ...
       a_range(opt_ind), chi2_a(opt_ind), P_values(opt_ind));

% Given the p-values of the direct method and pseudo-experimental method,
% data is not "fully" described by this optimal background combination.
% However, combination aMC1 + (1-a)MC2 describes data much better
% than the MC2 only, based on P-values.

