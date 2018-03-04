% Statistical Methods in Physics Spring 2014
% University of Helsinki
%
% Exercise 1
% Mikael Mieskolainen
clear; close all;

global h;
h = 2.4;

rho = @(w,l) (h*(w+l)) ./ (sqrt(l.^2 + w.^2)*sqrt(2*h.^2 + (w+l).^2));

Z = [];
i = 1; j = 1;
W = 0.1:0.25:100;
L = 0.1:0.25:100;

for w = L
    i = 1;
    for l = L
        Z(i,j) = rho(w,l);
        i = i + 1;
    end
    j = j + 1;
end

% Plot correlelation coefficient path
imagesc(W,L,Z); axis square;
xlabel('w (m)'); ylabel('l (m)'); title('\rho(\sigma_A_C, \sigma_A_W)');
colorbar;
