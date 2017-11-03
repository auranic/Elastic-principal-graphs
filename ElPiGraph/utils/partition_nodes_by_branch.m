function [node_partition, internal_flag, star_centers] = partition_nodes_by_branch(ElasticMatrix)
%% Decompose the graph into branches (stretches of edges connecting star centers)
% Input: representation of the primitive elastic graph as an
%   ElasticMatrix. The matrix will be transformed into an adjacency matrix,
%   so it can be simple an adjacency matrix as well
% Output:
%   node_partition - vector of size(ElasticMatrix,1) length with a
%      numbers indicating distinct branches. Star centers do not belong
%      to branches (they are documented separately in star_centers output)
%   internal_flag - vector of size(ElasticMatrix,1) length with 0/1
%      elements indicating if a node belongs to an internal branch (connecting
%      two stars (1) or to external branch, connecting a leaf and a star 
%      or another leaf (0)
%   star_centers - vector of structures, each having one field 'center'
%      indicating the node in the cencer, and one field 'branches' with a
%      vector of branch ids connected to the star. The order of the ids can
%      be further made data-dependent, using 'order_branches_in_stars'
%      function, to be used further in graph layout algorithms
%%
    % Convert ElasticMatrix into adjacency matrix and create copy
    A = (ElasticMatrix-diag(diag(ElasticMatrix)))>0;
    B = A;
    % Define size, preallocate arrays and calculate connectivities
    N = size(A,1);
    node_partition = zeros(N, 1);
    internal_flag = node_partition;
    Connectivities = sum(A);

    % Select terminal nodes and nodes with two lwaves
    term = Connectivities < 3;
    termInd = find(term);
    count = 1;
    while true
        % Search the first such node which is not coloured yet
        inds1 = termInd(node_partition(term) == 0);
        % If all such nodes are coloured then stop.
        if isempty(inds1) 
            break; 
        end;
        % seeding the branch
        ind = false(1, N);
        ind(inds1(1)) = true;
        node_partition(ind) = count;
        % extending the branch to the 'left'
        while true
            % Select all nodes, connected with previously selected
            inds = sum(A(ind, :), 1)>0;
            % Clear this node connection in A
            A(:, ind)=0;
            A(ind, :)=0;
            % Put current branch code to selected nodes
            node_partition(inds) = count;
            % Remove all exclude term
            ind = inds & term;
            if sum(ind) == 0
                break;
            end
        end
        % Check the type of branch. If there are two stars then it is
        % internal branch
        if nargout > 1
            if sum(Connectivities(node_partition == count) > 2) == 2
                internal_flag(inds) = 1;
            end
        end
        count = count+1;
    end

    % Remove branch marks from stars
    node_partition(~term)=0;
    internal_flag(~term)=0;

    if nargout < 3
        return;
    end
    
    % Find all stars with more than 2 leaves.
    stars = find(~term);
    if isempty(stars)
        star_centers = struct([]);
        return
    end
    
    for i=length(stars):-1:1
        inds = B(stars(i), :) & term;
        branches = node_partition(inds);
        sc = struct('center', stars(i), 'branches', branches);
        star_centers(i)=sc;
    end
end

function [A,node_partition] = extend_branch(A,node_partition,connectivities,k,branch_id)
    % if a star center or a leaf then stop extending
    node_partition(k) = branch_id;
    if connectivities(k)>2
        return;
    else
        inds = find(A(:,k)>0);
        A(:,k)=0;
        A(k,:)=0;
        for i=1:size(inds,1)
            [A,node_partition] = extend_branch(A,node_partition,connectivities,inds(i),branch_id);
        end
    end
end
