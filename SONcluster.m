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

%
% Make some simulated data.
%
seed = 12345;
stream = RandStream('mt19937ar','Seed',seed);
RandStream.setGlobalStream(stream);
X = 0.25*randn(n,p);
X(:,1:25) = X(:,1:25) + ones(2,25);

plot(X(1,:),X(2,:),'.')
%%
for i=1:p
    S(:,i) = X(:,i) - reshape(sum(L(i,:,:),2) - sum(L(:,i,:),1),n,1) ...
        + reshape(rho*(sum(V(i,:,:),2) - sum(V(:,i,:),1)),n,1); 
end

U = S*(1/nu)*(eye(p,p) - (2*rho/(1-4*rho*p))*ones(p,p));

for i=2:p
    for j=1:i-1
        V(i,j,:) = U(:,i) - U(:,j) + (1/rho)*reshape(L(i,j,:),n,1);
        V(j,i,:) = V(i,j,:);
    end
end

for i=2:p
    for j=1:i-1
        L(i,j,:) = reshape(L(i,j,:),n,1) + rho*(U(:,i) - U(:,j) - reshape(V(i,j,:),n,1));
        L(j,i,:) = L(i,j,:);
    end
end