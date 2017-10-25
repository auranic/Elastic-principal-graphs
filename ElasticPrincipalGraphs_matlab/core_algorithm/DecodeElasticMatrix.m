function [Edges,Lambdas,Mus] = DecodeElasticMatrix(ElasticMatrix)
%
% Converts ElasticMatrix into a set of edges and vectors of elasticities
% for edges and stars
%

Mus = diag(ElasticMatrix);
L = ElasticMatrix - diag(Mus);
[row,col] = find(L);
Edges = [row,col];

inds = find(Edges(:,1)<Edges(:,2));

Edges = Edges(inds,:);
Lambdas = L(Edges,1);
Lambdas = Lambdas(inds);

end
