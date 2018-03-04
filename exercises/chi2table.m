% Create chi2 table
%
clear;
close all

MAXNDF = 10;

for df = 1:MAXNDF
    C(df,:) = chi2inv([0.68 0.95],df);
end

disp([(1:MAXNDF)', C])
