function ElasticMatrix = Encode2ElasticMatrix(Edges, Lambdas, Mus)
% Converts array of edges and vectors of elasticities for edges and stars
% into ElasticMatrix 
%
%Inputs:
%   Edges is n-by-2 array with number of one edge's node in column 1 and
%       number of another node of the same edge in the column 2.
%   Lambdas is column vector of lambdas for edges specified in the Edges.
%   Mus is column vector of Mus specified for each node.
%

    %Define sizes
    NumberOfNodes = max(Edges(:));
    % Create ElasticMatrix
    ElasticMatrix = zeros(NumberOfNodes, NumberOfNodes);
    % Transform matrix of Edges into index of ElasticMatrix
    ind = sub2ind(size(ElasticMatrix), Edges(:,1), Edges(:,2));
    % Fill all edge elasticities
    ElasticMatrix(ind) = Lambdas;
    ElasticMatrix = ElasticMatrix + ElasticMatrix';
    % Add Mus to main diagonsl
    ElasticMatrix = ElasticMatrix+diag(Mus);
end