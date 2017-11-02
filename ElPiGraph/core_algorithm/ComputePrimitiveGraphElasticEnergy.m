function [ElasticEnergy, MSE, EP, RP] = ...
    ComputePrimitiveGraphElasticEnergy(NodePositions, ElasticMatrix,...
    dists, ~)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%   dists is n-by-1 vector. dists(i) is squared distance from data point
%       X(i,:) to node partition(i) or trimmed value.
%   BranchingFee is now unused fee for branching control (depicted as ~)
%
%Outputs
%   ElasticEnergy is total elastic energy 
%   MSE is mean square error of data approximation.
%   EP is edge potential 
%   RP is harmonicity potential 

    % Calculate MSE by usage dists
    MSE = sum(dists) / size(dists, 1);

    % Decompose ElasticMatrix
    Mu = diag(ElasticMatrix);
    Lambda = triu(ElasticMatrix,1);
    StarCenterIndices = find(Mu>0);
    
    % Calculate edge potential
    [row,col] = find(Lambda);
    dev = NodePositions(row,:)-NodePositions(col,:);
    l = Lambda(Lambda>0);
    EP = sum(l(:).*sum(dev.^2,2));

    % Form indeces
    indL = (Lambda + Lambda') > 0;

    % Calculate harmonicity potential
    RP = 0;
    for i=1:size(StarCenterIndices,1)
        leaves = indL(:,StarCenterIndices(i));
        K = sum(leaves);
        dev = NodePositions(StarCenterIndices(i),:)...
            - sum(NodePositions(leaves, :)) / K;
        RP = RP + Mu(StarCenterIndices(i)) * sum(dev .^ 2);
    end
    % Toatal energy is the sum.
    ElasticEnergy = MSE + EP + RP;
end
