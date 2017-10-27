function [NodeClusterCenters, NodeClusterRelativeSize] =...
    ComputeWeightedAverageM(X, partition, PointWeights, NumberOfNodes)
%ComputeWeightedAverage calculate NodeClusterCentres as weighted averages
%of points from matrix X.
%
%Inputs
%   X is n-by-m matrix of data points where each row corresponds to one
%       observation.
%   partition is n-by-1 (column) vector of node numbers. This vector
%       associate data points with Nodes.
%   PointWeights is n-by-1 (column) vector of point weights.
%   NumberOfNodes is number of nodes to calculate means.
%
%Important! if there is no point associated with node then coordinates of
%this node centroid are zero.
%
    % This recalcultion can be removed from this position
    X = bsxfun(@times, X, PointWeights);
    % Auxiliary calculations
    M = size(X,2);
    part = partition + 1;
    % Calculate total weights
    TotalWeight = sum(PointWeights);
    % Calculate weights for Relative size
    tmp = accumarray(part, PointWeights, [NumberOfNodes + 1, 1]);
    NodeClusterRelativeSize = tmp(2:end) / TotalWeight;
    
    NodeClusterCenters = zeros(NumberOfNodes + 1,size(X, 2));
    for k=1:M
        NodeClusterCenters(:, k) = accumarray(part,X(:, k), [NumberOfNodes+1, 1]) ./ tmp;
    end
    NodeClusterCenters = NodeClusterCenters(2:end,:);
end
