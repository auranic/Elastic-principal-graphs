function [graph] = makeGraph(NodePositions,Edges)

javaclasspath({'VDAOEngine.jar'});
graph = vdaoengine.analysis.grammars.Graph;

vdaoengine.analysis.grammars.Utils.setNodesForGraph(graph,NodePositions);
vdaoengine.analysis.grammars.Utils.setEdgesForGraph(graph,Edges);
graph.defineStarsFromPrimitiveGraphStructure();
graph.compileNodesInArrays();


