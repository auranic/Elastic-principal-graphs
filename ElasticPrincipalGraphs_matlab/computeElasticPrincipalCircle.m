function [NodePositions, Edges, ReportTable] = computeElasticPrincipalCircle(X,NumNodes,varargin)
%% computes principal circle (closed principal curve)
%

% we initialize the circle by placing 4 nodes in the plane of first two principal components

np = zeros(4,size(X,2));
[v,u,s] = pca(X);
mn = mean(X);
v1 = v(:,1)/norm(v(:,1));
v2 = v(:,2)/norm(v(:,2));
st1 = std(u(:,1));
st2 = std(u(:,1));
np(1,:) = mn-st1*v1'-st2*v2';
np(2,:) = mn-st1*v1'+st2*v2';
np(3,:) = mn+st1*v1'+st2*v2';
np(4,:) = mn+st1*v1'-st2*v2';
ed = [1,2;2,3;3,4;4,1];


[NodePositions,Edges,ReportTable] = computeElasticPrincipalGraph(X,NumNodes,...
    'GrowGrammars',[{'bisectedge'}],'ShrinkGrammars',[],...
    'InitGraph',struct('InitNodes',np,'InitEdges',ed),varargin{:});

end

