function [graph] = makeGraph(NodePositions,Edges)
%makeGraph make graph from node positions and edges
%   NodePositions is k-by-2 array of x and y coordinates of nodes.
%   Edges is k-by-2 matrix of integers. Edges(i,1) and Edges(i,2) specify
%       numbers of two vertex of i-th edge.
%Return object of vdaoengine.analysis.grammars.Graph class

    % Connection with Java to calculate
    javaclasspath({'.\core_algorithm_java\VDAOEngine.jar'});
    % Create graph
    graph = vdaoengine.analysis.grammars.Graph;
    % Initialise graph and extract it structure 
    vdaoengine.analysis.grammars.Utils.setNodesForGraph(graph,NodePositions);
    vdaoengine.analysis.grammars.Utils.setEdgesForGraph(graph,Edges);
    graph.defineStarsFromPrimitiveGraphStructure();
    graph.compileNodesInArrays();
end

