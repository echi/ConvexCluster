function [U,Vpath,Lpath] = SONcluster(X,lambdas,varargin)
%SONCLUSTER compute the SON cluster path on a set of points.
%
%   U = SONCLUSTER(X,lambdas) computes the cluster path for a matrix whose columns
%   are points to be clustered.
%
%   U = SONCLUSTER(X, lambdas, 'param', value, ...) specifies optional parameters and
%   values. Valid parameters and their default values are:
%      'tol' - Tolerance on relative change in the cluster centers {1e-8}
%      'maxiters' - Maximum number of iterations {1000}
%      'p' - which norm to use {2}
%      'gamma' - parameter that controls the local kernel weights {2}
%      'rho' - ADMM step size {10}
%      'printitn' - Print every n iterations; 0 for no printing {0}
%
% Author: Eric Chi



%% Set algorithm parameters from input or by using defaults.
params = inputParser;
params.addParamValue('tol',1e-8,@isscalar);
params.addParamValue('maxiters',1e3,@(x) isscalar(x) & x > 0);
params.addParamValue('norm',2,@(x) x==1 | x==2 | x==Inf);
params.addParamValue('gamma',2,@isscalar);
params.addParamValue('rho',10,@isscalar);
params.addParamValue('printitn',0,@isscalar);
params.parse(varargin{:});

%% Copy from params object
tol = params.Results.tol;
maxiters = params.Results.maxiters;
gamma = params.Results.gamma;
pnorm = params.Results.norm;
rho = params.Results.rho;
printItn = params.Results.printitn;

%% Obtain dimensions of the problem.
[d,N] = size(X);
q = N*(N-1)/2;
nLambda = numel(lambdas);
U = cell(nLambda,1);
Vpath = cell(nLambda,1);
Lpath = cell(nLambda,1);

%% Initialize variables
Q = zeros(q,N);
W = ones(q,1);
key = zeros(q,2);

i=1;
for t=1:N
    for k=t+1:N
        Q(i,t)=1;
        Q(i,k)=-1;
        key(i,:) = [t,k];
        W(i) = exp(-gamma*norm(X(:,t)-X(:,k),2)^2);
        i=i+1;
    end
end

if pnorm==1
    prox = @(z,tau) prox_L1(z,tau);
elseif pnorm==2
    prox = @(z,tau) prox_L2(z,tau);
elseif pnorm==Inf
    prox = @(z,tau) prox_Linf(z,tau);
end
    
V = zeros(d,q);
L = zeros(d,q);

nu = 1/(1 + N*rho);
u = X;
u0 = nu*X + (1-nu)*repmat(mean(X,2),1,N);

for ilambda=1:nLambda
    lambda = lambdas(ilambda);
    for iter=1:maxiters
%% Update U
        u_last = u;
        u = u0 + nu*(rho*V-L)*Q;
        if (norm(u-u_last)/(1+norm(u_last)) < tol)
            break;
        end
        UQ = u*Q';
%% Update V        
        Z = UQ + (1/rho)*L;
        V = prox(Z,(lambda/rho)*W');

%% Update L
        L = L + rho*(UQ-V);
        if (mod(iter,printItn)==0)
            fprintf(' Iter %4d: Objective = %.6e\n', ...
            iter, 0.5*norm(X-u)^2 + lambda*sum(norms(UQ,pnorm,2)));
        end        
    end
    U{ilambda} = u;
    Vpath{ilambda} = V;
    Lpath{ilambda} = L;
    fprintf('Lambda %d completed.\n',ilambda);
end