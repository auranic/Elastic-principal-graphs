function [Edges, Lambdas, Mus] = DecodeElasticMatrix(ElasticMatrix)
%
% Converts ElasticMatrix into a set of edges and vectors of elasticities
% for edges and stars
%

    Mus = diag(ElasticMatrix);
    L = ElasticMatrix - diag(Mus);
    [row, col] = find(triu(ElasticMatrix,1));
    Edges = [row,col];
    Lambdas = L(Edges,1);
end
