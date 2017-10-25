function [ElasticEnergy, MSE, EP, RP] = ...
    ComputePrimitiveGraphElasticEnergy(X, NodePositions, ElasticMatrix,...
    partition, dists, BranchingFee)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Computes elastic energy of primitive elastic graph 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%Inputs
%   X is n-by-m matrix of datapints with one data point per row. n is
%       number of data points and m is dimension of data space.
%   NodePositions is k-by-m matrix of embedded coordinates of graph nodes,
%       where k is number of nodes and m is dimension of data space.
%   ElasticMatrix is k-by-k matrix of nodes connectivity: 
%       ElsticMatrix(i,i) > 0 if node i is centre of star and zero otherwise
%       ElsticMatrix(i,j) > 0 if there is edge between nodes i and j. In
%       this case ElsticMatrix(i,j) is elasticity modulo of edge from i to j.
%   partition is n-by-1 vector. partition(i) is number of node which is
%       associated with data point X(i,:).
%   dists is n-by-1 vector. dists(i) is squared distance from data point
%       X(i,:) to node partition(i) or trimmed value.
%   BranchingFee is now unused fee for branching ?????
%
%Outputs
%   ElasticEnergy is total elastic energy (link???)
%   MSE is mean square error of data approximation.
%   EP is edge potential (link ???)
%   RP is harmonicity potential (link ???)

% Calculate MSE by usage dists
MSE = sum(dists)/size(X,1);

EP = 0;
RP = 0;


Mu = diag(ElasticMatrix);
Lambda = ElasticMatrix - diag(Mu);
StarCenterIndices = find(Mu>0);

[row,col] = find(Lambda);

for i=1:size(row,1)
    dev = NodePositions(row(i),:)-NodePositions(col(i),:);
    l = Lambda(row(i),col(i));
    EP = EP+l*(dev*dev');
end


for i=1:size(StarCenterIndices,1)
    leafs = find(Lambda(:,StarCenterIndices(i))>0);
    K = size(leafs,1);
    dev = NodePositions(StarCenterIndices(i),:)-1/K*sum(NodePositions(leafs,:));
    RP = RP+Mu(StarCenterIndices(i))*(dev*dev');
end

ElasticEnergy = MSE+EP+RP;

end

