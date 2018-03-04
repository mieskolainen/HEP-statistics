%% Task 2
% Statistical Methods
% Home Exam, Spring 2014
%
% Mikael Mieskolainen

% .avi files are saved animations of the boarding simulations

clear; close all;


%% i.) and ii.)

%{
N = 100; % The number of simulation runs

animation = false;
save_video = false;

% Simulation with outside-in boarding
bt1 = plainsimu(N, 'outside-in', animation, save_video);

% Simulation with random boarding
bt2 = plainsimu(N, 'random', animation, save_video);

save('./matlabdata/boardingtimes.mat', 'bt1', 'bt2');
%}
load('./matlabdata/boardingtimes.mat');

fprintf('Task 2. i.) Outside-in boarding time mean = %0.1f s, std = %0.1f s \n', ...
        mean(bt1), std(bt1));

fprintf('Task 2. ii.) Random boarding time mean = %0.1f s, std = %0.1f s\n', ...
        mean(bt2), std(bt2));

    