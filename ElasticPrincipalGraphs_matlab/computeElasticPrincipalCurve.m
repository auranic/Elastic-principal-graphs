function [NodePositions, Edges, ReportTable] = computeElasticPrincipalCurve(X,NumNodes,varargin)
%% computes principal curve
%
% The arguments are the same as in computeElasticPrincipalGraph
% Of note: if 'InitGraph' option is used, then this function will "grow"
% the initial graph structure by inserting nodes inside edges, thus
% smoothing the initial structure
%

[NodePositions,Edges,ReportTable] = computeElasticPrincipalGraph(X,NumNodes,'GrowGrammars',[{'bisectedge'}],'ShrinkGrammars',[],varargin{:});

end
