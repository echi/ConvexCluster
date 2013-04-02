%%
clear;

N = 50;
d = 2;

seed = 12345;
stream = RandStream('mt19937ar','Seed',seed);
RandStream.setGlobalStream(stream);
x = 0.05*randn(d,N);
x(:,1:25) = x(:,1:25) + ones(2,25);



Q = zeros(N*(N-1)/2,N);
key = zeros(N*(N-1)/2,2);

i=1;
for t=1:N
    for k=t+1:N
        Q(i,t)=1;
        Q(i,k)=-1;
        key(i,:) = [t,k];
        i=i+1;
    end
end

%%
close all;
figure;

%%
lambda_start = 0;
lambda_end = max(norms(x*Q',Inf,1))/N;
lambdas = lambda_start:0.00125:lambda_end;
nLambda = numel(lambdas);
loss = zeros(nLambda,1);

F(nLambda) = struct('cdata',[],'colormap',[]);
% Record the movie
plot(x(1,:),x(2,:),'b.'); hold on; grid;
F(1) = getframe;
for j = 2:nLambda
    lambda = lambdas(j);
    cvx_begin
    variable mu1(d,N)
    minimize(0.5*sum(sum((x-mu1).*(x-mu1))) ...
         +lambda*sum(norms(Q*mu1',1,2)))
    cvx_end
    plot(mu1(1,:),mu1(2,:),'r.');
    title(sprintf('lambda = %g',lambda));
    F(j) = getframe;
    loss(j) = 0.5*sum(norms(x - mu1,2,1).^2) + lambdas(j)*sum(norms(Q*mu1',2,d));
end

%%

rho = 10;
lambda =  lambdas(10);
eps_rel = 1e-13;
eps_abs = 1e-13;

q = size(Q,1);
V = zeros(d,q);
L = zeros(d,q);

nu = 1/(1 + N*rho);
tol = 1e-8;
max_iter = 2e4;
r = zeros(max_iter,1);
s = zeros(max_iter,1);
U = x;
U0 = nu*x + (1-nu)*repmat(mean(x,2),1,N);
for iter=1:max_iter
%% Update U
    U_last = U;
%    R = x + (rho*V-L)*Q;
    U = U0 + nu*(rho*V-L)*Q;
%    U = nu*R + (1-nu)*repmat(mean(R,2),1,N);
      if (norm(U-U_last)/(1+norm(U_last)) < tol)
          break;
      end
    UQ = U*Q';
%% Update V
    Z = (rho/lambda)*UQ + (1/lambda)*L;
    nz = norms(Z,2,1);
    z = (nz - 1)./nz;
    z(find(z < 0)) = 0;
    V = (lambda/rho)*bsxfun(@times,Z,z);

%% Update L
    L_last = L;
    L = L + rho*(UQ-V);
    
%% Primal residual
    r(iter) = norm(UQ - V,'fro');
%% Dual residual
    s(iter) = rho*norm((L - L_last)*Q,'fro');
    eps_pri = sqrt(d*N*(N-1)/2)*eps_abs + eps_rel*max(norm(UQ,'fro'),norm(V,'fro'));
    eps_dual = sqrt(d*N)*eps_abs + eps_rel*norm(L*Q,'fro');
%      if (r(iter) < eps_pri && s(iter) < eps_dual)
%          break;
%      end
end
r = r(1:iter);
s = s(1:iter);

%% Clustering assignment
C = cell(N,numel(V));
for i=1:numel(V)
    ix = find(norms(V{i},2,1) == 0);
    A = sparse(key(ix,1),key(ix,2),ones(numel(ix),1),N,N);
    visited = zeros(N,1);
    k = 1;

    for j=1:N
        if (visited(j) == 0)
            C{k,i} = find(bfs(A,j) >= 0);
            visited(C{k,i}) = 1;
            k = k + 1;
        end
    end
end

%% Compute cluster centroids
centroids = cell(N,numel(V));
for i=1:numel(V)
   for j=1:N
      if (~isempty(C{j,i}))
        centroids{j,i} = mean(x(:,C{j,i}),2);
      end
   end
end

%% Compute objective for difference lambda values
scores = zeros(numel(V),1);
den = zeros(numel(V),1);
xbar =  mean(x,2);
for i=1:numel(V)
   for j=1:N
      if (~isempty(C{j,i}))
        scores(i) = scores(i) + sum(norms(x(:,C{j,i}) - repmat(centroids{j,i},1,numel(C{j,i})),2,2).^2);
        den(i) = den(i) + numel(C{j,i})*norm(xbar-centroids{j,i},2)^2;
      end
   end
   scores(i) = scores(i)/den(i);
end

figure
plot(x(1,:),x(2,:),'b.'); hold on; grid on;
plot(U(1,:),U(2,:),'r.');

%%
F(nLambda+1) = struct('cdata',[],'colormap',[]);
% Record the movie
plot(x(1,:),x(2,:),'b.'); hold on; grid on;
F(1) = getframe;
for j = 1:nLambda
    plot(U{j}(1,:),U{j}(2,:),'r.');
    title(sprintf('lambda = %g',lambdas(j)));
    F(j+1) = getframe;
end
