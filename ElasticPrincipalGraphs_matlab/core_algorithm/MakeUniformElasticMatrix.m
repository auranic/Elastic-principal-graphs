function ElasticMatrix = MakeUniformElasticMatrix(Edges, Lambda, Mu)
% Creates an ElasticMatrix for a primitiveGraph defined by its vector of
%            edges
% The elasticities of the ElasticMatrix are equal to Lambda for all edges
% and to Mu for all stars.
%
% Inputs:
%   Edges is n-by-2 array with number of one edge's node in column 1 and
%       number of another node of the same edge in the column 2.
%   Lambda is edge elasticity.
%   Mu is star elasticity.
%
    %Define sizes
    NumberOfNodes = max(Edges(:));
    % Create ElasticMatrix
    ElasticMatrix = zeros(NumberOfNodes, NumberOfNodes);
    % Transform matrix of Edges into index of ElasticMatrix
    ind = sub2ind(size(ElasticMatrix), Edges(:,1), Edges(:,2));
    % Fill all edge elasticities
    ElasticMatrix(ind) = Lambda;
    ElasticMatrix = ElasticMatrix + ElasticMatrix';
    % Calculate connectivities
    Connect = sum(ElasticMatrix > 0);
    % Identify stars
    ind  = Connect > 1;
    % Generate vector for diagonal
    Mus = zeros(NumberOfNodes, 1);
    % Specify star elasticities
    Mus(ind) = Mu;
    % Add Mus to main diagonsl
    ElasticMatrix = ElasticMatrix + diag(Mus);
end

