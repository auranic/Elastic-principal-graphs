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



%[NodePositions,Edges,ReportTable] = computeElasticPrincipalGraph(X,NumNodes,...
%    'TrimmingRadius',TrimmingRadius,...
%    'InitGraph',struct('InitNodes',np,'InitEdges',ed),varargin{:});

end

