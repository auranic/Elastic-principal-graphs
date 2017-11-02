function [neighbour_nodes] = GetNeighbourhoodOnTheGraph(ElasticMatrix,nodes,radius)
% The function gives a set of nodes connected to the 'nodes' input by at most 
% 'radius' connections
% version 1.1

neighbour_nodes = nodes;
L = ElasticMatrix - diag(diag(ElasticMatrix));
L=L>0;

for i=1:radius
    addednodes = [];
    for j=1:size(neighbour_nodes,1)
        x = L(neighbour_nodes(j),:);
        nnodes = find(x>0);
        nnodes = nnodes';
        addednodes = cat(1,addednodes,nnodes);
    end
    list = cat(1,neighbour_nodes,addednodes);
    %% this is a fast implementation of unique for integers
    neighbour_nodes = find(accumarray(list(:)+1,1))-1;
end


end

