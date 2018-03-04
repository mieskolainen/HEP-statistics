% MLCG random number generator
% Multiplicative (m = 0) Linear Congruential Generator
%
% http://en.wikipedia.org/wiki/Linear_congruential_generator
%
% INPUT:    n  =  number of random numbers to be generated
%           s  =  seed (optional), something between
%                 0 <= s <= m, where here m = 2^31 - 1
%
% OUTPUT:   r  =  random number vector (n x 1)
%
% Mikael Mieskolainen

function r = MLCG(n,s)

% Fixed starting seed, 0 <= s <= m
if (nargin == 1)
    s = 2^31 / 10; % ad-hoc
end

% Parameters, see Wikipedia
m = 2^31 - 1; % the "modulus"
a = 48271;    % Multiplier
c = 0;        % the increment

% Output vector
r = zeros(n,1);

% Recurrence relation loop
for i = 1:n
    
    % Draw random number and normalize by dividing with m to [0,1]
    s = mod((a*s + c), m);
    r(i) = s / m;
end

end