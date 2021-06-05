function [NodeClusterCenters,NodeClusterRelativeSize] = ComputeWeightedAverage(X, partition, PointWeights, NumberOfNodes)
%ComputeWeightedAverage calculate NodeClausterCentres as weighted averages
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

    n = size(X,2);
    [rind, cind] = ndgrid(partition, 1:n);
    % Calculate sums.
    NodeClusterCenters = accumarray([rind(:), cind(:)],...
        reshape(bsxfun(@times,X,PointWeights),[],1),[NumberOfNodes,n]);
    % Calculate sum of weights
    NodeClusterRelativeSize = accumarray(partition, PointWeights);
    % Normalise NodeClusterCenters
    ind = NodeClusterRelativeSize > 0;
    NodeClusterCenters(ind,:) = bsxfun(@rdivide, NodeClusterCenters(ind,:),...
        NodeClusterRelativeSize(ind));
    % Normalise NodeClusterRelativeSize
    NodeClusterRelativeSize = NodeClusterRelativeSize...
        / sum(NodeClusterRelativeSize);
end