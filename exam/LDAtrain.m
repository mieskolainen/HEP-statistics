% Fisher's Linear Discriminant Classifier for the 2 class case
%
% NOTE! Function assumes that the priors are included already in the
% training set, i.e. the class fractions ~ priors
%
% INPUT:    X  =  Training data (d x N), where d is the data dimension
%           T  =  Class targets vector (N x 1), 0 = class A, 1 = class B
%
% OUTPUT:   w  =  Weight vector of the LDA
%           mu =  Cell array of class mean vectors
%           C  =  Cell array of class covariance matrices
%           
% Mikael Mieskolainen

function [w,mu,C] = LDAtrain(X,T)

% Class mean vectors
mu{1} = mean(X(:,T == 0), 2);
mu{2} = mean(X(:,T == 1), 2);

% Class covariance matrices
C{1} = cov(X(:,T == 0)');
C{2} = cov(X(:,T == 1)');

% Calculate Fisher LDA statistic, '\' is the inverse of C
w = (C{1} + C{2}) \ (mu{1} - mu{2});

end