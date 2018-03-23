function [partition, dists] = ...
    PartitionData(X, NodePositions, MaxBlockSize, SquaredX, TrimmingRadius)
%%%%%%%%%%%%%%%%%%%%%%%
%% Partition the data by proximity to graph nodes (same step as in K-means EM procedure)
%%%%%%%%%%%%%%%%%%%%%%%
%
%Inputs:
%   X is n-by-m matrix of datapoints with one data point per row. n is
%       number of data points and m is dimension of data space.
%   NodePositions is k-by-m matrix of embedded coordinates of graph nodes,
%       where k is number of nodes and m is dimension of data space.
%   MaxBlockSize integer number which defines maximal number of
%       simultaneously calculated distances. Maximal size of created matrix
%       is MaxBlockSize-by-k, where k is number of nodes.
%   SquaredX is n-by-1 vector of data vectors length: SquaredX = sum(X.^2,2); 
%   TrimmingRadius (optional) is squared trimming radius.
%
% Outputs
%   partition is n-by-1 vector. partition(i) is number of node which is
%       associated with data point X(i,:).
%   dists is n-by-1 vector. dists(i) is squared distance between node with
%       number partition(i) and data point X(i,:). 
%
    if nargin < 5
        TrimmingRadius = Inf;
    end
    if nargin < 4 || isempty(SquaredX)
        SquaredX = sum(X .^ 2, 2);
    end
    if nargin <3 || MaxBlockSize == 0
        MaxBlockSize = floor(10000000 / size(NodePositions, 1));
    end
    
    n = size(X, 1);
    partition = zeros(n, 1);
    dists = zeros(n, 1);
    %Calculate squared length of centroids
    cent = NodePositions';
    centrLength = sum(cent.^2);
    %Process partitioning without trimming
    for i = 1:MaxBlockSize:n
        % Define last element for calculation
        last = i + MaxBlockSize - 1;
        if last > n
            last = n;
        end
        % Prepare index
        ind = i:last;
        % Calculate distances
        d = bsxfun(@minus, centrLength,  2 * (X(ind,:) * cent));
        [dists(ind), partition(ind)] = min(d,[],2);
    end
    dists = dists + SquaredX;
    %Apply trimming
    ind = dists > TrimmingRadius;
    partition(ind) = 0;
    dists(ind) = TrimmingRadius;
end