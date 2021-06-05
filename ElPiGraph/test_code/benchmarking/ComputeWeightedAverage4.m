function [NodeClusterCenters,NodeClusterRelativeSize] = ComputeWeightedAverage4(X, partition, PointWeights, NumberOfNodes)
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
    NodeClusterCenters = zeros(NumberOfNodes,size(X,2));
    NodeClusterRelativeSize = zeros(NumberOfNodes,1);
    TotalWeight = sum(PointWeights);
    for i=1:NumberOfNodes
        inds = find(partition==i);
        if(size(inds,1)>0)
            tmp = sum(PointWeights(inds));
            NodeClusterCenters(i,:) = sum(bsxfun(@times,X(inds,:),PointWeights(inds)))/tmp;
            NodeClusterRelativeSize(i) = tmp/TotalWeight;
        end
    end
end