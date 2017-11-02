function [ElasticMatrix] = MakeUniformElasticMatrix(Edges,Lambda,Mu)
% Creates an ElasticMatrix for a primitiveGraph defined by its vector of
%            edges
% The elasticities of the ElasticMatrix are equal for all edgesa and Stars

NumberOfNodes = max(max(Edges));
NumberOfEdges = size(Edges,1);

ElasticMatrix = zeros(NumberOfNodes,NumberOfNodes);
AdjacencyMatrix = zeros(NumberOfNodes,NumberOfNodes);

for i=1:NumberOfEdges
    ElasticMatrix(Edges(i,1),Edges(i,2)) = Lambda;
    ElasticMatrix(Edges(i,2),Edges(i,1)) = Lambda;
    AdjacencyMatrix(Edges(i,1),Edges(i,2)) = 1;
    AdjacencyMatrix(Edges(i,2),Edges(i,1)) = 1;
end

Connectivities = sum(AdjacencyMatrix);
MuVector = Mu*(Connectivities>1);

ElasticMatrix = ElasticMatrix+diag(MuVector);
end

