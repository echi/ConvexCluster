function y = prox_L2(x,tau)
%PROX_L computes the proximal mapping associated with tau*norm(x,2)
%
%   y = PROX_L2(x,tau) comptues the proximal mapping associated with
%   2-norm.
%
% Author: Eric Chi
z = 1 - tau./norms(x,2,1);
z(z < 0) = 0;
y = bsxfun(@times,x,z);
%y = max(1-tau/norm(x,2),0)*x;
