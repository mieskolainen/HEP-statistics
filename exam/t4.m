%% Task 4
% Statistical Methods
% Home Exam, Spring 2014
%
% Mikael Mieskolainen

clear; close all;


%% i.)

% The spesific heat capacity (C) of a substance is the amount of heat energy
% (E) required to change its temperature by a given amount per unit mass
% (m), i.e. C = E / (m deltaT), where deltaT is the change in temperature.
m = 1000; % mass of liquid

% Electric heater produces heat as E = Pt, where power P = U x I is
% constant as a function of time with
U = 12;   % Volts
I = 10;   % Amperes

% Time points and their symmetric  +- errors/uncertainties
t = [352 701 1048 1398 1751 2099 2446 2805];
t_err = [5 9 9 9 9 15 15 15];

% Temperature change and their symmetric +- errors/uncertainties
deltaT = [10.0 19.7 30.2 40.4 49.9 60.5 70.4 80.0];
deltaT_err = [0.1 0.2 0.2 0.2 0.3 0.3 0.4 0.4];

% Let's fit a linear regression line y = Cx using least squares here
% E = C*m*deltaT <=>
% U*I*t = C*m*deltaT <=>
% t = C*(m/(U*I))*deltaT  ~ t = y and x = (m/(U*I))*deltaT
t_func = @(deltaT, C, m, U, I) C*(m/(U*I))*deltaT;

% Thus, uncertainties are neglected in x (in deltaT ~ temperature measur.)
% because it's an implicit assumption in the ordinary least squares. Errors
% exist only in y-variable, i.e. in time (t).

% Model matrix (column vector)
X = m/(U*I)*deltaT';

% Response vector (y) and its error, and error matrix
y = t';
y_err = t_err';
V = diag(y_err .^2);

% LS estimate of the parameter, and it's variance using pseudoinverse
% formula
C_hat = (X' * inv(V) * X) \ X' * inv(V) * y;
C_var = inv(X' * inv(V) * X);

fprintf('Task 4. i.) The heat capacity C = %0.3f +- %0.3f (Joule/g*K) \n', ...
        C_hat, sqrt(C_var));

% Chi^2 in matrix formalism
chi2 = (y - X*C_hat)' / V * (y - X*C_hat);

% Degrees of freedom = # measurements - # parameters to fit
DOF = length(y) - 1;

% P-value
P_value = 1 - chi2cdf(chi2, DOF);
fprintf('Task 4. i.) chi^2: %0.3f (P-value: %0.3f) \n', chi2, P_value); 

% Plot out our fit and measurements
figure;
deltaT_range = 0:0.01:90;
plot(deltaT_range, t_func(deltaT_range, C_hat, m, U, I), 'k-'); hold on;
errorbar(deltaT, t, t_err, 'r.');
xlabel('\DeltaT (Celsius)'); ylabel('t (s)');
legend(sprintf('LS-fit, \\chi^2 = %0.3f', chi2), 'Data');


%% ii.)

% Let's check if the variance estimate obtained in i.) is at the RCF-bound.
% Well, we know already that it is because it was a linear least squares,
% which is efficient and we used a closed-form estimator for the variance
% of the parameter estimate. Let's do it anyway,
% by differentiating negative log-likelihood in the RCF-formula twice with
% respect to parameter C evaluated at C = C_hat, i.e.
% -d^2lnL/dCdC, where ln L = -chi^2 / 2
%                          = -1/2 x sum_i ( t_i - Cx_i)^2 / sigma_t_i^2
% 
% (G. Cowan Statistical Data Analysis)
%
% and taking the inverse of the result gives us:
C_var_rcf = 1 / sum(X.^2 ./ t_err.^2');

fprintf('Task 4. ii.) Variance in i.) %0.6f, RCF-bound %0.6f \n', ...
        C_var, C_var_rcf);

% So, our estimator in i.) is at the RCF-bound.


%% iii.)

% Now taking into account errors in both x and y,
% no closed form estimator in general,
% find the minimum using a brute-force search

% Range of brute force scan, hand tuned
C_range = linspace(C_hat - 0.02, C_hat + 0.02, 100);

%{
% Random initializations, very sensitive to this, must be large!
N = 1e5;

% Chi2 values
chi2c = zeros(length(C_range), 1);

% Thetas (parameter C estimates)
thetas = zeros(length(C_range), 1);

% x-vector and its error (NOTE, here include the m,U and I factors, which
% is a must because our equation is basically y = Cx)
x = m/(U*I)*deltaT';
x_err = m/(U*I)*deltaT_err';

% Create random samples by variating x-vector
X_new = zeros(length(deltaT),N);
for k = 1:N
    X_new(:,k) = x + x_err.*randn(size(x));
end

i = 1;
for C1 = C_range
    chi2_min = inf;

    for k = 1:N
        x_new = X_new(:,k);

        % Chi2 of "Total Least Squares" type of cost function,
        % errors in both x and y
        chi2  = sum( ( y - x_new*C1 ).^2 ./ y_err.^2 ) ...
              + sum( ( x_new - x ).^2 ./ x_err.^2 );
    
        if (chi2 < chi2_min)
            chi2_min = chi2;
        end
    end
    chi2c(i) = chi2_min;
    thetas(i) = C1;
    i = i + 1;
end

% Save variables
save('./matlabdata/chi2loop.mat', 'chi2c', 'thetas');

%}

% Load variables
load('./matlabdata/chi2loop');

% Find the minimum among all parameter pairs
chi2_min = inf;
theta_min = 0;
min_ind = 0;

for i = 1:size(chi2c,1)
    if (chi2c(i) < chi2_min)
        chi2_min = chi2c(i);
        theta_min = thetas(i);
        min_ind = i;
    end
end

% Find the errors from the chi^2 plot "~ Graphical Method"
% A one sigma deviation is at the point(s) given by:
% chi2_min + deviation (see e.g. G. Cowan, Statistical Data Analysis)
C1_err = [0,0];

deviation = 1; % chi^2 + 1.0

% Find lower side error
for i = min_ind(1):-1:1
   if (chi2c(i) >= (chi2_min + deviation))       
       C1_err(1) = theta_min - thetas(i);
       break;
   end
end
% Find upper side error
for i = min_ind(1):size(chi2c,1)
   if (chi2c(i) >= (chi2_min + deviation))
       C1_err(2) = thetas(i) - theta_min;
       break;
   end
end

% Plot chi^2 distribution and chi2_min + deviation error bounds
figure;
plot(C_range, chi2c); hold on;
plot(theta_min - [C1_err(1) C1_err(1)], [min(chi2c) max(chi2c)], 'k--');
plot([theta_min theta_min], [min(chi2c) max(chi2c)], 'k-');
plot(theta_min + [C1_err(2) C1_err(2)], [min(chi2c) max(chi2c)], 'k--');
plot(C_range, ones(size(C_range))*(chi2_min + deviation), 'r--');
axis square; axis tight;
xlabel('C');
ylabel('\chi^2 distribution');
title(sprintf('Parameter C lower and upper \\chi^2 + %0.1f errors (dashed black)', deviation));

% P-value
P_value = 1 - chi2cdf(chi2_min, DOF);

fprintf('Task 4. iii.) Min(chi^2): %0.3f (P-value: %0.3f) \n', ...
        chi2_min, P_value);
fprintf('Task 4. iii.) Parameter C = %0.3f with errors (-%0.3f, +%0.3f) \n', ...
        theta_min, C1_err(1), C1_err(2));

% The change in chi^2 and P-value is towards slightly better fit, i.e.
% smaller chi^2 and larger P-value, which is reasonable because we take
% into account errors both in x and y. Though, not much of change in the
% parameter C value estimate.


%% iv.)

% The chemist's error estimates are a little large in time (t) variable,
% I think they could be much smaller in reality, because those are order
% of a wall clock errors...

% The temperature measurement errors seem quite ok, given e.g. some
% semiconductor based measurement sensor.

% Given the quite large P-value ~ 0.34 (or in iii.) ~ 0.64 ),
% results are compatible with the theory, i.e. no need to neglect null
% hypothesis (H_0 ~ our function matches the data). Also, because
% P-value is reasonable (not extremely high), the errors are probably
% estimated quite ok.

