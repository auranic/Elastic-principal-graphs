function neighbour_nodes =...
    GetNeighbourhoodOnTheGraph(ElasticMatrix, nodes, radius)
% The function gives a set of nodes connected to the 'nodes' input by at
% most 'radius' edges on graph (nodes with shortest path to 'nodes' which
% is not longer than 'radius'. 
% version 1.1
%
% Inputs:
%   ElasticMatrix - matrix with elasticity coefficients (lambdas for edges,
%       and star elasticities along the diagonal
%   nodes is list of nodes whose neighbours have to be found
%   radius is radius of neighbourhood.
%
% Output is list of nodes in neighbourhood of radius 'radius' for nodes
% 'nodes'.

    % Create vector for index
    neighbours = false(1, size(ElasticMatrix, 1));
    neighbours(nodes) = true;
    % Form sceleton of Elastic matrix
    L = ElasticMatrix > 0;
    % Add neighbour of radius r
    for r = 1:radius
        neighbours = neighbours | (sum(L(neighbours, :), 1)>0);
    end
    neighbour_nodes = find(neighbours);
end

