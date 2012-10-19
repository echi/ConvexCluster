%
% SON clustering
%
% Date: October 10, 2012
% Author: Eric Chi

rho = 2;
p = 50;
n = 2;
nu = 1 - 2*p*rho;

L = zeros(p,p,n);
V = zeros(p,p,n);
X = zeros(n,p);
S = zeros(n,p);

for i=1:p
    S(:,i) = X(:,i) - reshape(sum(L(i,:,:),2) - sum(L(:,i,:),1),n,1) ...
        + reshape(rho*(sum(V(i,:,:),2) - sum(V(:,i,:),1)),n,1); 
end
