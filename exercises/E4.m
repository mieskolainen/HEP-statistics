% Statistical Methods in Physics Spring 2014
% University of Helsinki
%
% Exercise 4
% Mikael Mieskolainen

addpath ./kolmofunc % Kolmogorov distribution functions

clear;
close all;


%% TASK 1

% Read in data for class A
fid = fopen('./data/sample_cata.txt');
XA = fscanf(fid, '%f %f', [2 inf]);
fclose(fid);

% Read in data for class B
fid = fopen('./data/sample_catb.txt');
XB = fscanf(fid, '%f %f', [2 inf]);
fclose(fid);

% Class targets, 0 == class A, 1 == class B
T = [zeros(size(XA,2),1); ones(size(XB,2),1)];

% Train the LDA, get weight vector w (the c in exercise sheet notation)
X = [XA XB];
[w,mu,C] = LDAtrain(X, T);

% The weight vector
w


%% i.)

% Plot data points in 2D
figure;
plot(XA(1,:), XA(2,:), 'ro'); hold on;
plot(XB(1,:), XB(2,:), 'xb');
xlabel('x_1'); ylabel('x_2');
legend('Class A', 'Class B');


%% ii.)

figure;
set(gcf,'visible','off');

% Calculate the projection (test statistic of LDA) for each event
t = X'*w;

% Plot distributions
t_bins = -21:0.075:-3;
[nA, ~] = hist(t(T == 0), t_bins);
pdfA = nA / sum(nA); % Normalized to pdf
bar(t_bins, pdfA, 'r'); hold on;

[nB, ~] = hist(t(T == 1), t_bins);
pdfB = nB / sum(nB); % Normalized to pdf
bar(t_bins, pdfB, 'b');

xlabel('$t = \langle \mathbf{x},\mathbf{w} \rangle$', ...
       'interpreter', 'latex');
ylabel('pdf(t)', 'interpreter', 'latex');


%% iii.)

% Calculate cumulative distributions
CDF_A = cumsum(pdfA);
CDF_B = cumsum(pdfB);

% Calculate the t_cut necessary to get 95% rejection from class A
% (find the bin which is the closest)
[~, ind] = min(abs(CDF_A - 0.05));
t_cut = t_bins(ind);
fprintf('T1 iii.) Condition t_{cut} = %0.1f, gives 95%% rejection of class A \n', ...
        t_cut);

% Calculate the efficiency of class B with this t_cut
fprintf('T1 iii.) Efficiency of class B with t_{cut} = %0.1f is %0.2f%% \n', ...
        t_cut, CDF_B(ind));

% Plot t_cut in the previous figure
set(gcf,'visible','on');
plot([t_cut t_cut], [0 0.06], 'k');

legend('Class A', 'Class B', 't_{cut}'); axis tight;
    
% Lets plot the events with this selection rule t < t_cut
%{
figure;
X_cutA = X(:, t < t_cut & T == 0);
X_cutB = X(:, t < t_cut & T == 1);

plot(X_cutA(1,:), X_cutA(2,:), 'ro'); hold on;
plot(X_cutB(1,:), X_cutB(2,:), 'xb');
xlabel('x_1'); ylabel('x_2');
legend('Class A', 'Class B');
%}

% Plot data points in 2D
figure;
plot(XA(1,:), XA(2,:), 'ro'); hold on;
plot(XB(1,:), XB(2,:), 'xb');
xlabel('x_1'); ylabel('x_2');

hold on;

% Affine decision boundary plot
x = linspace(min(X(1,:)), max(X(1,:)), 2);
% <x,w> = t_cut <=> x2 = (t_cut - x1w1)w2
d = [x; (t_cut - x*w(1))/w(2)];
plot(d(1,:), d(2,:), '-k'); hold on;

% Plot orthogonal projection plane (line)
theta = pi/2; % rotate 90 deg
R = [cos(theta) -sin(theta); 
     sin(theta) cos(theta)];
d_ = R*d;
xy_off = [1, -0.5]; % Arbitrary for visual purposes
plot(d_(1,:) + xy_off(1), d_(2,:) + xy_off(2), '--k');

legend('Class A', 'Class B', ...
       'Decision boundary: <x,w> = t_{cut}', ...
       'Fisher projection plane (t-space)', 'Location', 'SouthEast');
axis equal;


%% TASK 2


%% i.) and ii.)

% The recoil angles 
data1 = [80, 44, 56, 65, 33, 52, 42, 104];
data2 = [18, 108, 70, 135, 68, 32, 30, 38, 122, 49, 91];

% PDFs and CDFs
bins = 0:1:180;
[n1,~] = hist(data1, bins); 
[n2,~] = hist(data2, bins);
pdf1 = n1 / sum(n1);
pdf2 = n2 / sum(n2);
cdf1 = cumsum(pdf1);
cdf2 = cumsum(pdf2);

figure;
stairs(bins, cdf1, 'r--'); hold on;
stairs(bins, cdf2, 'b'); axis tight;
xlabel('\theta (deg)');
ylabel('CDF(\theta)');
legend('Data1', 'Data2', 'Location', 'NorthWest');

% Kolmogorov-Smirnov test
KS_value = max(abs(cdf1 - cdf2));
fprintf('T2 i.) KS-value %0.4f \n', KS_value);


% The null hypothesis H0 ~ both data come from the same distribution
% The alternative hypothesis Ha ~ data come from different distributions

% Significance level:  \alpha = 0.05
% Critical value: D_alpha
% Critical region: Reject H0 if D > D_alpha

% You reject the null hypothesis 
% that the two samples were drawn from the same
% distribution if the p-value is less than your
% significance level \alpha.

N1 = length(data1);
N2 = length(data2);

P = 1 - kolmcdf(KS_value * sqrt((N1*N2)/(N1 + N2)));
fprintf('T2 ii.) P-value %0.4f \n', P);

% P_value more than \alpha = 0.05, we don't reject H_0 so
% they are probably coming from the same distribution

% Test the MATLAB function
%[h,p,ks2test] = kstest2(data1, data2, 0.05);


%% iii.) and iv.)

% Null hypothesis H0 ~ theoretical distribution and data compatible
% Alternative hypothesis H_a ~ data and theoretical incompatible

% Theta range (180 ... 0 deg), cos(theta) = -1 ... 1
N_bins = 250;
costheta = linspace(-1,1,N_bins);

% Theoretical distribution dN/dcos(theta) ~ (up to some constant)
% The effect of constant is normalized away in the CDF [0,1]
gamma = 0.1;
% Actually if gamma = 1, data fits well with prediction, try it :D

%f = @(costheta) 1 + gamma*costheta;

% Integral formula analytically: N = int (1 + gamma*costheta) dcos(theta)
% =>
g = @(costheta) costheta + gamma*0.5*costheta.^2;

% Calculate integrals to get the theoretical CDF
Z = zeros(size(costheta));

% NOTE that these integrals must be uniformly spaced on cos(theta)
% domain [-1,...,1], because we do our KS-comparison there, if one
% would integrate here with uniform stepping in the angle theta domain,
% the integral would be wrongly nonlinearly "warped" at cos(theta) domain.
% (Calculus 101 - change of variables changes also the integration domain!)
for i = 2:length(costheta)
    interval = [costheta(1)  costheta(i)];
    %Z(i) = trapz(interval, f(interval));   % Numerical integr.
    Z(i) = g(interval(2)) - g(interval(1)); % Analytical integr.
end

% Normalize the CDF
theory_CDF = Z / Z(end);

% Take cosine of the data for the comparison: theta -> cos(theta)
D1 = cos(deg2rad(data1));
D2 = cos(deg2rad(data2));

% Calculate normalized CDF for data1 + data2
% We give the bin edges, so use histc instead of hist
[n,~] = histc([D1 D2], costheta);
PDF12 = n / sum(n);
CDF12 = cumsum(PDF12);

% Plot theory and data
figure;
stairs(costheta, theory_CDF, 'b--'); hold on;
stairs(costheta, CDF12, 'r'); 
stairs(costheta, abs(theory_CDF - CDF12), 'k-.');
axis tight;
xlabel('cos(\theta)');
ylabel('CDF(cos(\theta))');
legend('Theory', 'Data1 + Data2', '| CDF_{theory} - CDF_{data}|', ...
       'Location', 'NorthWest');

KS_value2 = max(abs(theory_CDF - CDF12));
fprintf('T2 iii.) KS-value %0.4f \n', KS_value2);

% iv.)

% P-value is
P = 1 - kolmcdf(KS_value2 * sqrt(N1 + N2));
fprintf('T2 iv.) P-value %0.4f \n', P);

% P-value less than chosen significance level alpha = 0.05, we reject H_0,
% data and theoretical distribution not compatible!

