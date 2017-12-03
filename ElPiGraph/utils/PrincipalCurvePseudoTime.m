function [coordinate,traversal,X_projected] = PrincipalCurvePseudoTime(X,nodes,edges,terminal)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    %ElasticMatrix = Encode2ElasticMatrix(edges,1,1);
    %A = (ElasticMatrix-diag(diag(ElasticMatrix)))>0;
    N = max(max(edges));    
    A = zeros(N,N);
    for i=1:size(edges,1)
        A(edges(i,1),edges(i,2)) = 1;
        A(edges(i,2),edges(i,1)) = 1;
    end
    Connectivities = sum(A);
    traversal = zeros(1,N);

    ind = terminal;
    count = 1;
    
    traversal(count) = ind;
    
    while true
        ind1 = find(sum(A(ind, :), 1)>0);
        A(:, ind)=0;
        A(ind, :)=0;
        if isempty(ind1) 
            break; 
        end;
        count=count+1;
        traversal(count) = ind1;
        ind = ind1;
    end

    SquaredX = sum(X .^ 2, 2);
    partition = PartitionData(X, nodes, 100000, SquaredX);
    [MSE, X_projected, EdgeIndices, ProjectionValues] = project_point_onto_graph(X, nodes, edges, partition);
    
    coordinate = zeros(size(X,1),1);
    for i=1:size(X,1)
        e = edges(EdgeIndices(i),:);
        v = ProjectionValues(i);
        e1 = find(traversal==e(1,1));
        e2 = find(traversal==e(1,2));
        if(e1<e2) 
            coordinate(i) = e1+v;
        else
            coordinate(i) = e2+(1-v);
        end
    end
    
    coordinate = (coordinate-min(coordinate))/(max(coordinate)-min(coordinate));

end

