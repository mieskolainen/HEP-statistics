% Statistical Methods in Physics Spring 2014
% University of Helsinki
%
% Exercise 7
% Mikael Mieskolainen

clear; close all;


%% Task 1

% T = 2 pi sqrt (L/g)
% where T is the pendulum period, L is the length and g is the
% gravitational acceleration

% Fit a model y = theta*x
% theta ~ 2 pi / sqrt(g)
% x ~ sqrt(L)
% y ~ T

% Model matrix (in meters)
X = [sqrt(0.85) sqrt(0.75) sqrt(0.70) sqrt(0.60) sqrt(0.50)]';

% Measurement vector y
y = [1.850 1.737 1.678 1.554 1.419]';

% Measurement error matrix
sigma2 = [0.012 0.012 0.011 0.010 0.010].^2;
V = diag(sigma2);


%% i.)

% LS estimate
theta_hat = (X' * inv(V) * X) \ X' * inv(V) * y;
theta_cov = inv(X' * inv(V) * X);

% Transform of variable
g_hat = ((2*pi)/theta_hat)^2;

% Error propagation to get error of g
f = @(theta, sigma2) sqrt( ( (-8*pi^2)/(theta^3) )^2 * sigma2 );
g_err = f(theta_hat, theta_cov);

fprintf('Task 1: i.) g_hat = %0.3f +- %0.3f m/s^2\n', g_hat, g_err);

% Chi^2
chi2 = (y - X*theta_hat)' / V * (y - X*theta_hat);

% Degrees of freedom = # measurements - # parameters to fit
DOF = length(y) - 1;

% P-value
P_value = 1 - chi2cdf(chi2, DOF);

fprintf('Task 1: i.) chi^2: %0.4f (P-value: %0.2f) \n', chi2, P_value); 

%% ii.)

% Calculate error using RFC as theta = theta_hat
% Formula below obtained by differentiating twice

theta_var_rfc = 1 / sum(X.^2 ./ sigma2');
g_err_rfc = f(theta_hat, theta_var_rfc);

fprintf('Task 1: ii.) Error of g using RFC is %0.3f \n', g_err_rfc);

% Same as using "normal/matrix equations", i.e. Maximum Likelihood solution
% is the same as Least Squares in the case of Gaussian measurement noise!


%% iii.)

% P-value is 1, i.e. fit is really super. There might be some
% alteration of results, especially if these are hand made "student
% measurements". Given errors/uncertainties on y are also
% especially small ... :D


%% Task 2

warning off; % Remove singularity warnings

% Model
% y = a1 x - a2 / x

% y - Measurements
y = [-4.02 -2.74 -1.15 1.49 6.87]';
sigma2_y = ([0.50 0.25 0.08 0.09 1.90].^2)';

% Error matrix for least squares (errors only in y)
V = diag(sigma2_y);

% x - Measurements
x = [22000, 22930, 23880, 25130, 26390]';
sigma2_x = ([440 470 500 530 540].^2)';

% Model matrix
X = [x -1./x];

%% i.)

% LS fit (x is accurate i.e. no errors,
% which is an implicit asumption always in the ordinary LS method!)
theta_hat = (X' / V * X) \ X' / V * y;
theta_cov = inv(X' / V * X);

% Chi^2
chi2 = (y - X*theta_hat)' / V * (y - X*theta_hat);

% Degrees of freedom = # measurements - # parameters to fit
DOF = length(y) - 2;

% P-value
P_value = 1 - chi2cdf(chi2, DOF);

fprintf('Task 2: i.) chi^2: %0.3f (P-value: %0.3f) \n at (a1, a2) = (%0.6f +- %0.6f, %0.3f +- %0.3f) \n', ...
 chi2, P_value, theta_hat(1), sqrt(theta_cov(1,1)), ...
 theta_hat(2), sqrt(theta_cov(2,2)) );

% System parameters
C = 0.02e-6; % Farads
w0 = 1;      % rad/s

% Calculate resistance (Ohms) and its error via error prop.
R = 1/(w0*C*theta_hat(2));
R_err = sqrt( (-1/(w0*C*theta_hat(2)^2))^2*theta_cov(2,2) );
fprintf('Task 2: i.) R = %0.2f +- %0.2f (Ohm)\n', R, R_err);

% Calculate inductance (Henry) and its error via error prop.
% ~ (df/da)^2*sigma_a^2 + (df/db)^2*sigma_b^2 + 2 df/da df/db * cov_ab 
% L = (theta_hat(1)*R) / w0;
L = theta_hat(1) / (C*w0^2*theta_hat(2));
L_err = sqrt( (1/(C*w0^2*theta_hat(2)))^2*theta_cov(1,1) + ...
              (theta_hat(1)/(C*w0^2*theta_hat(2)^2))^2*theta_cov(2,2) + ...
              (2*theta_hat(1))/(C^2*w0^4*theta_hat(2)^3)*theta_cov(1,2));

fprintf('Task 2: i.) L = %0.3f +- %0.3f (H)\n', L, L_err);


%% ii.)

% Find minimum using brute-force search

% Range of brute force scan, hand tuned
a1_range = linspace(theta_hat(1)/2, 1.5*theta_hat(1), 100);
a2_range = linspace(theta_hat(2)/2, 1.5*theta_hat(2), 100);

%{
tic;

% Random initializations
N = 12000;

% Chi2 contour
chi2c = zeros(length(a1_range), length(a2_range));

% Thetas
thetas = cell(length(a1_range), length(a2_range));

% Create random samples by variating x
X_new = zeros(length(x),N);
for k = 1:N
    X_new(:,k) = x + sqrt(sigma2_x).*randn(size(x));
end

i = 1;
for a1 = a1_range
    j = 1;
    for a2 = a2_range
        chi2_min = inf;

        for k = 1:N
            
            x_new = X_new(:,k);

            % Chi2 of "total least squares" type of cost function,
            % errors in both x and y
            chi2  = sum( ( y - [x_new -1./x_new]*[a1; a2] ).^2 ./ sigma2_y ) ...
                  + sum( ( x_new - x ).^2 ./ sigma2_x );
            
            if (chi2 < chi2_min)
                chi2_min = chi2;
            end
        end
        
        chi2c(i,j) = chi2_min;
        thetas{i,j} = [a1; a2];
        j = j + 1;
    end
    i = i + 1;
end

toc;

% Save variables
save('chi2loop', 'chi2c', 'thetas');

%}

% Load variables
load('chi2loop');



%% Find the minimum among all parameter pairs

chi2_min = inf;
theta_min = [0;0];
min_ind = [0;0];

for i = 1:size(chi2c,1)
    for j = 1:size(chi2c,2)
        if (chi2c(i,j) < chi2_min)
            chi2_min = chi2c(i,j);
            theta_min = thetas{i,j};
            min_ind = [i;j];
        end
    end
end

figure;
imagesc(a1_range, a2_range, log10(chi2c)); colorbar;
axis square; axis tight;
xlabel('\alpha_1'); ylabel('\alpha_2');
title('log_{10}(\chi^2) distribution');


%% Find the errors from the chi^2 contour
% A one sigma deviation from the point(s):
% chi2_min + deviation. However, see the uncertainty coverage of
% multidimensional parameter estimation in G. Cowan, Statistical Data
% Analysis.

deviation = 1.0; % The basic

a1_err = [0,0];
a2_err = [0,0];

% Find lower side error in a1
for i = min_ind(1):-1:1
   if (chi2c(i,min_ind(2)) >= (chi2_min + deviation))       
       a1_err(1) = theta_min(1) - thetas{i,min_ind(2)}(1);
       break;
   end
end
% Find upper side error in a1
for i = min_ind(1):size(chi2c,1)
   if (chi2c(i,min_ind(2)) >= (chi2_min + deviation))
       a1_err(2) = thetas{i,min_ind(2)}(1) - theta_min(1);
       break;
   end
end
% Find lower side error in a2
for j = min_ind(2):-1:1
   if (chi2c(min_ind(1),j) >= (chi2_min + deviation))
       a2_err(1) = theta_min(2) - thetas{min_ind(1),j}(2);
       break;
   end
end
% Find upper side error in a2
for j = min_ind(2):size(chi2c,2)
   if (chi2c(min_ind(1),j) >= (chi2_min + deviation))
       a2_err(2) = thetas{min_ind(1),j}(2) - theta_min(2);
       break;
   end
end

% Average upper and lower errors and get square => variance
a1_var = mean(a1_err)^2;
a2_var = mean(a2_err)^2;

% Put into parameter covariance matrix (neglect covariance estimation)
theta_cov = [a1_var 0;
             0 a2_var];

% Degrees of freedom = # measurements - # parameters to fit
DOF = 5 - 2;

% P-value
P_value = 1 - chi2cdf(chi2_min, DOF);

fprintf('Task 2: ii.) Minimum chi^2: %0.3f (P-value: %0.3f) at \n (a1, a2) = (%0.6f +- %0.6f, %0.3f +- %0.3f) \n', ...
 chi2_min, P_value, theta_min(1), sqrt(a1_var), theta_min(2), sqrt(a2_var));


%% Parameter estimates for R and L


% Calculate resistance (Ohms) and its error via error prop.
R = 1/(w0*C*theta_min(2));
R_err = sqrt( (-1/(w0*C*theta_min(2)^2))^2*theta_cov(2,2) );
fprintf('Task 2: ii.) R = %0.2f +- %0.2f (Ohm)\n', R, R_err);

% Calculate inductance (Henry) and its error via error prop.
% (neglect covariance here)
% ~ (df/da)^2*sigma_a^2 + (df/db)^2*sigma_b^2 + 2 df/da df/db * cov_ab 
% L = (theta_hat(1)*R) / w0;
L = theta_min(1) / (C*w0^2*theta_min(2));
L_err = sqrt( (1/(C*w0^2*theta_min(2)))^2*theta_cov(1,1) + ...
              (theta_min(1)/(C*w0^2*theta_min(2)^2))^2*theta_cov(2,2));

fprintf('Task 2: ii.) L = %0.3f +- %0.3f (H)\n', L, L_err);

% Some changes in comparison with part i.), not drastical,
% the value of R bit lower, P-value larger => indicating better fit,
% and also uncertainty (error) estimates are smaller for L and R

