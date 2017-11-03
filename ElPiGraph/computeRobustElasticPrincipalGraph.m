function [NodePositions, Edges, ReportTable] = computeRobustElasticPrincipalGraph(X,NumNodes,TrimmingRadius,varargin)
%% computes robust principal graph (with trimming)
%

% we initialize the graph by placing 2 nodes close to each other, in the
% densest part of the data distribution

np = zeros(2,size(X,2));
ed = [1,2];

NumberOfSamples = floor(size(X,1)/10);
if NumberOfSamples>1000
    NumberOfSamples=1000;
end

sampling = randperm(size(X,1),NumberOfSamples);

[dense_point_index, neighbours] = find_highest_local_density_point(X, sampling);

np(1,:) = X(dense_point_index,:);
np(2,:) = X(neighbours(1),:);

NodePositions = np;
Edges = ed;

    [NodePositions,Edges,ReportTable] = computeElasticPrincipalGraph(X,NumNodes,...
    'TrimmingRadius',TrimmingRadius,...
    'InitGraph',struct('InitNodes',np,'InitEdges',ed),varargin{:});

end

