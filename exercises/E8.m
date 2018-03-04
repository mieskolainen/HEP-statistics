% Statistical Methods in Physics Spring 2014
% University of Helsinki
%
% Exercise 8
% Mikael Mieskolainen

clear; close all;


%% Task 2.

% Height (m)
z = [0 6 12 18]*1e-6;

% The number of clusters
n = [1880 940 530 305];

global dRho;
dRho = 63;         % Density difference kg/m^3
global r;
r = 0.52e-6;       % Radius of the sphere (m)
global g;
g = 9.80;          % m/s^2
global T;
T = 293;           % Temperature (K)

% Physical formula
v = @(z,v0,k) v0*exp(-4*pi*r^3*dRho*g*z/(3*k*T));


%% i.) Binned Maximum Likelihood for the parameter k

v0 = n(1); % Assume v0 = n (z = 0)

k_range = 0.8e-23:0.0001e-23:1.60e-23;
lnL = zeros(size(k_range));
ind = 1;

for k = k_range
    for i = 1:length(n)
        lnL(ind) = lnL(ind) + n(i)*log(v(z(i),v0,k)) - v(z(i),v0,k);
    end
    ind = ind + 1;
end

[max_lnL_value, max_lnL_ind] = max(lnL);
k_hat = k_range(max_lnL_ind);

% Uncertainty on k using: log(L) - 1/2
% However,
% see more about the multidimensional uncertainty region coverage in e.g.
% Cowan's book - not so simple as in 1D.

% Lower error
error_low = 0;
error_low_ind = 0;
for i = max_lnL_ind:-1:1
   if (lnL(i) <= max_lnL_value - 0.5) % LogL - 1/2
       error_low = k - k_range(i);
       error_low_ind = i;
       break;
   end
end

% Upper error
error_up = 0;
error_up_ind = 0;
for i = max_lnL_ind:length(lnL)
   if (lnL(i) <= max_lnL_value - 0.5) % LogL - 1/2
       error_up = k_range(i) - k;
       error_up_ind = i;
       break;
   end
end

plot(k_range, lnL); hold on;
plot(ones(2,1)*k_range(error_low_ind), [min(lnL) max(lnL)], 'r--');
plot(ones(2,1)*k_hat, [min(lnL) max(lnL)], 'k');
plot(ones(2,1)*k_range(error_up_ind), [min(lnL) max(lnL)], 'r--');
axis tight;
xlabel('Boltzman constant k (J/K)');
ylabel('ln L(k)');
title('Binned ML estimation of the parameter k');

% Take the error as average of lower and upper errors
k_err = mean([error_low error_up]);

fprintf('Task 2: i.) k = %0.2f +- %0.2f x 10^(-23) J/K \n', ...
        k_hat*1e23, k_err*1e23);

    
%% ii.) Binned Maximum Likelihood for the parameters v0 and k

k_range = 0.8e-23:0.0005e-23:1.60e-23;
v0_range = 1600:0.5:2200;

lnL = zeros(length(k_range), length(v0_range));

ind_i = 1;
for k = k_range
    ind_j = 1;
    for v0 = v0_range
        for i = 1:length(n)
            lnL(ind_i, ind_j) = lnL(ind_i, ind_j) ...
                              + n(i)*log(v(z(i),v0,k)) - v(z(i),v0,k);
        end
        ind_j = ind_j + 1;
    end
    ind_i = ind_i + 1;
end

% Find out the maximum values, i.e. our estimates
max_lnL = -inf;
max_inds = [0,0];
for i = 1:size(lnL,1)
    for j = 1:size(lnL,2)
        if (lnL(i,j) > max_lnL)
           max_lnL = lnL(i,j);
           max_inds = [i,j];
        end
    end
end
k_hat = k_range(max_inds(1));
v0_hat = v0_range(max_inds(2));

% Plot lnL surface
figure;
surface(v0_range, k_range, lnL, 'EdgeColor', 'none'); colorbar;
axis tight;
xlabel('v_0');
ylabel('k (J/K)');
title('ln L(v_0,k)');

% Plot estimate and 1 sigma confidence ellipsoid
figure;
plot(v0_hat, k_hat, 'ko'); hold on;
mask = lnL.*(lnL >= max(lnL(:)) - 0.5); % LogL - 1/2
contour(v0_range, k_range, mask, 'EdgeColor', 'black');
axis tight;
xlabel('v_0');
ylabel('k (J/K)');
title('1 sigma confidence ellipsoid');


% Find out the boundaries of the confidence ellipsoid by marginalizing
v0_mask = sum(mask,1) & ones(size(v0_range));
k_mask = sum(mask,2) & ones(size(k_range'));

% k errors
k_error_low = k_hat - min(k_range(k_mask == 1));
k_error_up = max(k_range(k_mask == 1)) - k_hat;

% v0 errors
v0_error_low = v0_hat - min(v0_range(v0_mask == 1));
v0_error_up = max(v0_range(v0_mask == 1)) - v0_hat;

% Take the error as the average of lower and upper errors
k_err = mean([k_error_low k_error_up]);
v0_err = mean([v0_error_low v0_error_up]);

fprintf('Task 2: ii.) k = %0.2f +- %0.2f x 10^(-23) J/K, v0 = %0.0f +- %0.0f \n', ...
        k_hat*1e23, k_err*1e23, v0_hat, v0_err);

%% iii.)

R = 8.314; % J/mol*K

% Avogadros constant
N_A = R / k_hat;

% The error of it using error propagation: sigma_N_A = R*sigma_k / k^2
N_A_err = R*k_err / k_hat^2;

fprintf('Task 2: iii.) N_A = %0.2f +- %0.2f x 10^(-23) 1/mol \n', ...
        N_A*1e-23, N_A_err*1e-23);

% Larger value obtained for N_A than what is known (because k is too small) 

