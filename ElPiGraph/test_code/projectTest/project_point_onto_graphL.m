function [MSE, X_projected, EdgeIndices, ProjectionValues] =...
    project_point_onto_graphL(X, NodePositions, Edges, partition)
%% This function calculates piece-wise linear projection of a dataset onto
%% the graph defined by NodePositions and Edges
% Input arguments:
%   X - is the n-by-m data matrix. Each row corresponds to one data point.
%   NodePositions - is k-by-m matrix of positions of the original graph
%       nodes in the same space as X.
%   Edges is k-by-2 matrix of integers. Edges(i,1) and Edges(i,2) specify
%       numbers of two vertex of i-th edge.
%   partition is n-by-1 vector. partition(i) is number of node which is
%       associated with data point X(i,:).
% 
% Outputs:
%   MSE is mean squared error (distance) from X to projection X onto graph
%   X_projected is matrix of the coordinates of data point projections
%   EdgeIndices is n-by-1 with number of edge which contains projection
%   ProjectionValues is n-by-1 real vector with values between 0 and 1
%       indicating where between Edge(index,1) and Edge(index,2) the
%       projection was done.  
%%

    %Preallocate arrays for result
    n = size(X, 1);
    X_projected = zeros(size(X));
    ProjectionValues = zeros(n,1);
    EdgeIndices = zeros(n,1);

    % Calculate Direction vectors for all edges
    vec = NodePositions(Edges(:, 2), :) - NodePositions(Edges(:, 1),:);
    vecSqL = sum(vec .^ 2, 2);
    
    % Auxiliary arrays
    best = Inf(n, 1);
    dis = Inf(n, 1);
    proj = zeros(n, 1);
    X_projected1 = zeros(size(X));
    
    % Go over the edges
    for ed = 1:size(Edges, 1)
        % Get index of points linked to this edge
        indP = partition == Edges(ed, 1) | partition == Edges(ed, 2);
        % Calculate projection
        u = bsxfun(@minus, X(indP, :), NodePositions(Edges(ed, 1), :))...
            * vec(ed, :)' / vecSqL(ed);
        u(u < 0) = 0;
        u(u > 1) = 1;
        % Find coordinates of projection
        xp = bsxfun(@plus, NodePositions(Edges(ed, 1),:), ...
            bsxfun(@times, vec(ed, :), u));
        % Calculate distance to projection
        xd = sum((X(indP, :) - xp) .^ 2, 2);
        % Put data to greater arrays
        dis(indP) = xd;
        proj(indP) = u;
        X_projected1(indP, :) = xp;
        % Select elements where calculated projection is the best
        indP = dis < best;
        best(indP) = dis(indP);
        X_projected(indP, :) = X_projected1(indP, :);
        ProjectionValues(indP) = proj(indP);
        EdgeIndices(indP) = ed;
    end
    MSE = mean(best);
end