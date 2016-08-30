function [LabelColorMap, graph, TaxonMap] = drawPieChartsMetroMap(NodePositions,Edges,data,labels,varargin)

if size(labels,1)>size(labels,2)
    display('ERROR: Labels argument must be row-vector');
else

nodesMM = computeMetroMapLayout(NodePositions,Edges);
graph = makeGraph(NodePositions,Edges);
TaxonMap = getTaxonMap(graph,data);

[LabelColorMap] = drawPieChartsProc(nodesMM, TaxonMap, labels, varargin{:});

end;

end