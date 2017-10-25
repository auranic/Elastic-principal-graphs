function [ElasticMatrix] = Encode2ElasticMatrix(Edges,Lambdas,Mus)
%
% Converts set of edges and vectors of elasticities
% for edges and stars into ElasticMatrix
%

NumberOfNodes = max(max(Edges));
NumberOfEdges = size(Edges,1);

ElasticMatrix = zeros(NumberOfNodes,NumberOfNodes);

for i=1:NumberOfEdges
    ElasticMatrix(Edges(i,1),Edges(i,2)) = Lambdas(i);
    ElasticMatrix(Edges(i,2),Edges(i,1)) = Lambdas(i);
end

ElasticMatrix = ElasticMatrix+diag(Mus);

end