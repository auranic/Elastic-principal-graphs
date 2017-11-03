function [X_projected, MSE, EdgeIndices, ProjectionValues] =...
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
%   X_projected is matrix of the coordinates of data point projections
%   MSE is mean squared error (distance) from X to projection X onto graph
%   EdgeIndices is n-by-1 with number of edge which contains projection
%   ProjectionValues is n-by-1 real vector with values between 0 and 1
%   indicating where between Edge(index,1) and Edge(index,2) the projection
%   was done. 
%%
    %Preallocate arrays for result
    N = size(X,1);
    X_projected = zeros(size(X));
    ProjectionValues = zeros(N,1);
    Distances_squared = zeros(N,1);

    % Calculate Direction vectors for all edges
    vec = NodePositions(Edges(:, 2), :)-NodePositions(Edges(:, 1),:);
    vec = bsxfun(@rdivide, vec, sum(vec .^ 2,2));
    

    for i=1:N
        k = partition(i);
        inds  = find((Edges(:,1)==k)|(Edges(:,2)==k));
        ds = zeros(size(inds,1),1);
        xp = zeros(size(inds,1),size(X,2));
        pv = zeros(size(inds,1),1);
        for j=1:size(inds,1)
            [xp(j,:),pv(j),ds(j)] = project_point_onto_edge(X(i,:),NodePositions,Edges(inds(j),:));
        end
        [mind,mini] = min(ds);
        X_projected(i,:) = xp(mini,:);
        ProjectionValues(i) = pv(mini);
        Distances_squared(i) = ds(mini);
    end

end

function [x_projected, projection_value, distance_squared] =...
    project_point_onto_edge(x, NodePositions, Edge)
    % Calculate 
vec = NodePositions(Edge(2),:)-NodePositions(Edge(1),:);
u = ((x-NodePositions(Edge(1),:))*vec')/(vec*vec');
if u<0
    x_projected = NodePositions(Edge(1),:);
elseif u>1
    x_projected = NodePositions(Edge(2),:);
else
    x_projected = NodePositions(Edge(1),:)+u*vec;
end

projection_value = u;
distance_squared = (x_projected-x)*(x_projected-x)';

end



