function [MSE, X_projected, EdgeIndices, ProjectionValues] =...
    project_point_onto_graph(X, NodePositions, Edges, partition)
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
    [n, m] = size(X);
    X_projected = zeros(size(X));
    ProjectionValues = zeros(n,1);
    EdgeIndices = zeros(n,1);
    MSE = 0;

    % Calculate Direction vectors for all edges
    vec = NodePositions(Edges(:, 2), :)-NodePositions(Edges(:, 1),:);
    vecSqL = sum(vec .^ 2, 2);
    
    for i=1:n
        k = partition(i);
        inds  = find((Edges(:,1)==k)|(Edges(:,2)==k));
        r = length(inds);
        ds = zeros(r, 1);
        xp = zeros(r, m);
        pv = zeros(r, 1);
        for j=1:size(inds,1)
            [xp(j,:),pv(j),ds(j)] = ...
                project_point_onto_edge2(X(i,:), NodePositions,...
                Edges(inds(j),:), vec(inds(j),:), vecSqL(inds(j)));
        end
        [mini, mInd] = min(ds);
        X_projected(i,:) = xp(mInd,:);
        ProjectionValues(i) = pv(mInd);
        EdgeIndices(i) = inds(mInd);
        MSE = MSE + mini;
    end
    MSE = MSE / n;
end

function [x_projected, u, distance_squared] =...
    project_point_onto_edge2(x, NodePositions, Edge, vec, vecL)
    % Calculate 
    u = ((x-NodePositions(Edge(1),:))*vec')/vecL;
    if u<0
        u = 0;
    elseif u>1
        u = 1;
    end
    x_projected = NodePositions(Edge(1),:) + u * vec;

    distance_squared = (x_projected-x)*(x_projected-x)';

end