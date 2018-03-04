% Statistical Methods in Physics Spring 2014
% University of Helsinki
%
% Exercise 6
% Mikael Mieskolainen

clear;
close all;

% Read in data
fid = fopen('./data/ml_sample_1.txt');
D20 = fscanf(fid, '%f', [1 inf])';
fclose(fid);

% Read in data
fid = fopen('./data/ml_sample_2.txt');
D100 = fscanf(fid, '%f', [1 inf])';
fclose(fid);


%% TASK 2

% PDF of data, t in [0,1]
pdf = @(f,t) abs(sin(2*pi*f*t));

% Frequency range uniformly in 0.1 ... 10
f = logspace(-1, 1, 100000);


%% i.) and iii.)

data = {D20, D100};
label = {'20', '100'};

% Estimated values
f_est = zeros(2,1);

% Estimated errors
f_err = zeros(2,2);

for k = 1:2

% Log likelihood as a function of f, prod(pdf)-> sum(log(pdf))
figure;

lnL = zeros(size(f));
for i = 1:length(f)
    lnL(i) = sum( log(pdf(f(i), data{k})) );
end
subplot(1,2,1);
semilogx(f, lnL);
set(gcf,'visible','off');
xlabel('f (Hz)');
ylabel(sprintf('ln L_{%s}(f)', label{k}));
axis square;

% Find the maximum ln L and its f value
[~, ind] = max(lnL);

fprintf('i / iii.) lnL_max,%s = %0.3f and f_%s = %0.3f (Hz) \n', label{k}, ...
        lnL(ind), label{k}, f(ind));

% Save value
f_est(k) = f(ind);


%% ii.) and iv.)

% Find the uncertainty using "graphical method"

target = lnL(ind) - 0.5;

% Lower side parabola scan, i.e. find where ln L_max - 1/2
lower_ind = 0;
diff_min = inf;
for i = ind:-1:1
    diff = abs(target - lnL(i));
    if (diff < diff_min)
        diff_min = diff;
        lower_ind = i;
    end
end

% Upper side parabola scan, i.e. find where ln L_max - 1/2
upper_ind = 0;
diff_min = inf;
for i = ind:length(lnL)
    diff = abs(target - lnL(i));
    if (diff < diff_min)
        diff_min = diff;
        upper_ind = i;
    end
end

% Now plot that
subplot(1,2,2);
semilogx(f(lower_ind:upper_ind), lnL(lower_ind:upper_ind));
set(gcf,'visible','on');
xlabel('f (Hz)');
ylabel(sprintf('ln L_{%s}(f)', label{k}));
axis square;

% The uncertainty corresponding to one standard deviation
% [f - lower_value, f + upper_value]
fprintf('ii / iv.) 1-sigma intervals: [f_%s - %0.3f, f_%s + %0.3f]\n', ...
        label{k}, f(ind) - abs(f(lower_ind)), label{k}, f(upper_ind) - f(ind));

% Save errors
f_err(k,:) = [f(ind) - abs(f(lower_ind)), f(upper_ind) - f(ind)];


end

%% v.)

% Taking the uncertainties

% Based on sample sizes, uncertainty should be ~ sqrt(20) / sqrt(100)
std_20 = f_est(1) / sqrt(20);
std_100 = f_est(2) / sqrt(100);
sqrtN_ratio = std_20 / std_100;

% Take the mean of lower and upper errors (they are really close)
graph_ratio = mean(f_err(1,:)) / mean(f_err(2,:));

fprintf('v.) Graphical errors ratio: %0.3f, sqrt(N) ratio: %0.3f \n', ...
        graph_ratio, sqrtN_ratio);

% Not exactly the same, but on the order of magnitude compatible

