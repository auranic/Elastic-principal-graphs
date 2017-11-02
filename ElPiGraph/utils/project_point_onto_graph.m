function [X_projected,MSE,EdgeIndices,ProjectionValues] = project_point_onto_graph(X,NodePositions,Edges,partition)
%% this function calculates piece-wise linear projection of a dataset onto
%% the graph defined by NodePositions and Edges
% Input arguments:
% X - is the dataset
% NodePositions, Edges - definition of graph embedment
% partition - integer vector of the length equal to size(X,1) defining
% the closest node index in the graph
% 
% Outputs:
% X_projected: projected dataset (same space as X)
% MSE - mean squared error (distance) from X to graph
% EdgeIndices - integer vector of length size(X,1) with indications of
% on which edges the projection was done
% ProjectionValues - real vector of length size(X,1) with values between 0
% and 1 indicating where between Edge(index,1) and Edge(index,2) the
% projection was done.

%%

%[X_projected,ProjectionValues, Distances_squared] = project_point_onto_edge(X,NodePositions,Edges);

X_projected = zeros(size(X,1),size(X,2));
ProjectionValues = zeros(size(X,1),1);
Distances_squared = zeros(size(X,1),1);



for i=1:size(X,1)
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

function [x_projected,projection_value, distance_squared] = project_point_onto_edge(x,NodePositions,Edge)
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



