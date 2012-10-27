function y = prox_L1(x,tau)
%PROX_L1 computes the proximal mapping associated with tau*norm(x,1)
%
%   y = PROX_L1(x,tau) comptues the proximal mapping associated with
%   1-norm.
%
% Author: Eric Chi
[m,n] = size(x);
tau = repmat(tau,m,1); % Assumes tau is a row vector of length n. Check it!
y = (abs(x) - tau);
y(y < 0) = 0;
y = sign(x).*y;