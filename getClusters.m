function [C, centroids, cluster_sizes] = getClusters(X,V)
%GETCLUSTERS takes the output of SONcluster and extracts the cluster
%assignments.
%
%   [C, centroids, cluster_sizes] = getClusters(X,V)
%
% Author: Eric Chi

N = size(X,2);
q = N*(N-1)/2;
nLambda = numel(V);
Q = zeros(q,N);
key = zeros(q,2);

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
C = cell(N,nLambda);
for i=1:nLambda
    ix = find(norms(V{i},2,1) == 0);
    A = sparse([key(ix,1);key(ix,2)],[key(ix,2);key(ix,1)],ones(2*numel(ix),1),N,N);
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

cluster_sizes = cellfun(@numel,C);
cluster_sizes = cluster_sizes(cluster_sizes > 0);

%% Compute cluster centroids
centroids = cell(N,nLambda);
for i=1:nLambda
   for j=1:N
      if (~isempty(C{j,i}))
        centroids{j,i} = mean(X(:,C{j,i}),2);
      end
   end
end
