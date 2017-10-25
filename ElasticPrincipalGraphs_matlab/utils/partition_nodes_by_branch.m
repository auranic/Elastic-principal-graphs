function [node_partition, internal_flag, star_centers] = partition_nodes_by_branch(ElasticMatrix)
%% Decompose the graph into branches (stretches of edges connecting star centers)
%   Input: representation of the primitive elastic graph as an
%   ElasticMatrix. The matrix will be transformed into an adjacency matrix,
%   so it can be simple an adjacency matrix as well
%   Output:
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

A = (ElasticMatrix-diag(diag(ElasticMatrix)))>0;
node_partition = zeros(size(A,1),1);
Connectivities = sum(A);

count = 1;
internal_flag = zeros(size(A,1),1);

while 1    
    inds = find((Connectivities==2)|(Connectivities==1));
    inds1 = [];
    l = 1;
    for j=1:size(inds,2)
        if node_partition(inds(j))==0 inds1(l)=inds(j); l=l+1; end;
    end
    if size(inds1,1)==0 break; end;
    % seeding the branch
    k = inds1(1);
    node_partition(k)=count;
    % extending the branch to the 'left'
    [A,node_partition] = extend_branch(A,node_partition,Connectivities,k,count);
    
    number_of_branches = max(node_partition);
    for i=1:number_of_branches
        inds = find(node_partition==i);
        number_of_stars_in_branch = sum(Connectivities(inds)>2);
        if number_of_stars_in_branch==2
            internal_flag(inds) = 1;
        end
    end
    count = count+1;
end
stars = find(Connectivities>2);
node_partition(stars)=0;
internal_flag(stars)=0;

if(size(stars,2)==0)
    star_centers = struct([]);
end

A = (ElasticMatrix-diag(diag(ElasticMatrix)))>0;
for i=1:size(stars,2)
    inds = find(A(stars(i),:)>0);
    branches = [];
    k=1;
    for j=1:size(inds,2) 
        if node_partition(inds(j))>0
            branches(k)=node_partition(inds(j)); k=k+1;
        end
    end
    sc = struct('center',stars(i),'branches',branches);
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
