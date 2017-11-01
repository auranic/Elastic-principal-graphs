function [NodePositionArray, ElasticMatrices, NodeIndicesArray] =...
    GraphGrammarOperation(NodePositions, ElasticMatrix, X, type, SqX)
%
% This is the core function for application of graph grammar approach for
% constructing primitive elastic principal graphs
% 
% The function takes a definition of the principal graph embeddment and
% applies a graph grammar operation of type type
% 
% Input arguments:
% NodePositions - position of nodes of the primitive elastic graph
% ElasticMatrix - matrix with elasticity coefficients (lambdas for edges,
% and star elasticities along the diagonal
% X - is dataset which is to be approximated by the graph
% type - one of the operation types:
%   'addnode2node'
%   'removenode'
%   'bisectedge'
%   'shrinkedge'
% SqX is column vector of squared lengths of data vectors
%
% Outputs:
% NodePositionArray - 3d array with dimensions
% [node_number, NodePosition, graph_number], represents all generated node
% configurations
% ElasticMatrices - 3d array with dimensions 
% [node_number, node_number, graph_number], 
% represents all generated elasticity matrices 
%
%% in the version 1.1 with a possibility of local search, each operation reports NodeIndices
%% which specifies how the nodes in the newly generated graph are related to the nodes
%% in the initial graph, and contains zeros for new nodes.

    switch type
        case 'addnode2node'
            [NodePositionArray, ElasticMatrices, NodeIndicesArray] = AddNode2Node(NodePositions,ElasticMatrix,X,SqX);
        case 'removenode'
            [NodePositionArray, ElasticMatrices, NodeIndicesArray] = RemoveNode(NodePositions,ElasticMatrix);
        case 'bisectedge'
            [NodePositionArray, ElasticMatrices, NodeIndicesArray] = BisectEdge(NodePositions,ElasticMatrix);
        case 'shrinkedge'
            [NodePositionArray, ElasticMatrices, NodeIndicesArray] = ShrinkEdge(NodePositions,ElasticMatrix);
        otherwise
            display(sprintf('ERROR: operation %s is not defined',type));
    end
end


function [NodePositionArray, ElasticMatrices, NodeIndicesArray]...
    = AddNode2Node(NodePositions, ElasticMatrix, X, SqX)
%
% This grammar operation adds a node to each graph node
% The positions of the node is chosen as a linear extrapolation for a leaf
% node (in this case the elasticity of a newborn star is chosed as in
% BisectEdge operation),
% or
% as the data point giving the minimum local MSE for a star (without any optimization).

    NNodes = size(NodePositions,1);
    NumberOfGraphs = NNodes;
    NodePositionArray = zeros(NNodes + 1, size(NodePositions, 2), NumberOfGraphs);
    ElasticMatrices = zeros(NNodes + 1, NNodes + 1, NumberOfGraphs);
    Mus = diag(ElasticMatrix);
    L = ElasticMatrix - diag(Mus);
    Connectivities = sum(L > 0);

    k=1;
    partition = PartitionData(X, NodePositions, 100000, SqX, Inf);
        
for i=1:NNodes
     
    meanLambda = mean(L(i,L(i,:)>0));

    if(Connectivities(i)==1)
        ineighbour = find(L(i,:)>0);
        NewNodePosition = 2 * NodePositions(i, :) - NodePositions(ineighbour, :);
        [np,em,inds] = f_add_nonconnected_node(NodePositions, ElasticMatrix, NewNodePosition);
        nn = size(np, 1);
        [em] = f_add_edge(em,i,nn,ElasticMatrix(i,ineighbour));
        em(i,i) = ElasticMatrix(ineighbour,ineighbour);
        NodePositionArray(:,:,k) = np(:,:);
        ElasticMatrices(:,:,k) = em(:,:);
        NodeIndicesArray(:,k) = inds;
        k=k+1;        
    else
        
        [nplocal,emlocal,inds] = f_get_star(NodePositions,ElasticMatrix,i);        
        %indlocal = [];
        %for l=1:size(inds,2)
        %   indlocal = [indlocal; find(partition==inds(l))];
        %end
         indlocal=find(partition==i);
         xlocal = X(indlocal,:);

        if size(xlocal,1)==0
            % empty star
            NodeNewPosition = mean(nplocal);
        else
            
        minMSE = realmax;
        m = -1;
        
        % mean point of the central cluster - seems to work the best
        NodeNewPosition = mean(xlocal);
        
        end
        
        [np,em,inds] = f_add_nonconnected_node(NodePositions,ElasticMatrix,NodeNewPosition);
        nn = size(np,1);
        [em] = f_add_edge(em,i,nn,meanLambda);
        NodePositionArray(:,:,k) = np(:,:);
        ElasticMatrices(:,:,k) = em(:,:);
        NodeIndicesArray(:,k) = inds;
        k=k+1;
    end
end

end


function [NodePositionArray, ElasticMatrices, NodeIndicesArray] = BisectEdge(NodePositions,ElasticMatrix)
%
% This grammar operation inserts a node inside the middle of each edge
% The elasticity of the edges do not change
% The elasticity of the newborn star is chosen as 
% mean over the neighbour stars if the edge connects two star centers
% or 
% the one of the single neigbour star if this is a dangling edge
% or 
% if one starts from a single edge, the star elasticities should be on
% one of two elements in the diagoal of the ElasticMatrix 

[Edges,Lambdas,Mus] = DecodeElasticMatrix(ElasticMatrix);
NumberOfGraphs = size(Edges,1);
NNodes = size(NodePositions,1);

NodePositionArray = zeros(NNodes+1,size(NodePositions,2),NumberOfGraphs);
ElasticMatrices = zeros(NNodes+1,NNodes+1,NumberOfGraphs);

k=1;
for i=1:size(Edges,1)
   NewNodePosition = (NodePositions(Edges(i,1),:)+NodePositions(Edges(i,2),:))/2;
   [np,em, inds] = f_add_nonconnected_node(NodePositions,ElasticMatrix,NewNodePosition);
   nn = size(np,1);
   lambda = ElasticMatrix(Edges(i,1),Edges(i,2));
   em = f_removeedge(em,Edges(i,1),Edges(i,2));
   em = f_add_edge(em,Edges(i,1),nn,lambda);
   em = f_add_edge(em,Edges(i,2),nn,lambda);
   mu1 = ElasticMatrix(Edges(i,1),Edges(i,1));
   mu2 = ElasticMatrix(Edges(i,2),Edges(i,2));
   if mu1>0&mu2>0
       em(nn,nn) = (mu1+mu2)/2;
   else
       em(nn,nn) = max(mu1,mu2);
   end
   NodePositionArray(:,:,k) = np(:,:);
   ElasticMatrices(:,:,k) = em(:,:);
   NodeIndicesArray(:,k) = inds;
   k=k+1;
end

end


function [NodePositionArray, ElasticMatrices, NodeIndicesArray] = RemoveNode(NodePositions,ElasticMatrix)
%
% This grammar operation removes a leaf node (connectivity==1)
%

Mus = diag(ElasticMatrix);
L = ElasticMatrix - diag(Mus);
Connectivities = sum(L>0);
NNodes = size(NodePositions,1);

NumberOfGraphs = sum(Connectivities==1);
NodePositionArray = zeros(NNodes-1,size(NodePositions,2),NumberOfGraphs);
ElasticMatrices = zeros(NNodes-1,NNodes-1,NumberOfGraphs);

k=1;
for i=1:size(Connectivities,2)
    if(Connectivities(i)==1)
        [np,em,inds] = f_remove_node(NodePositions,ElasticMatrix,i);
        NodePositionArray(:,:,k) = np(:,:);
        ElasticMatrices(:,:,k) = em(:,:);
        NodeIndicesArray(:,k) = inds;
        k=k+1;
    end
end

end

function [NodePositionArray, ElasticMatrices, NodeIndicesArray] = ShrinkEdge(NodePositions,ElasticMatrix)
%
% This grammar operation removes an edge from the graph
% If this is an edge connecting a leaf node then it is equivalent to
% RemoveNode. So we remove only internal edges.
% If this is an edge connecting two stars then their leaves are merged,
% and the star is placed in the middle of the shrinked edge.
% The elasticity of the new formed star is the average of two star
% elasticities.
%

Mus = diag(ElasticMatrix);
L = ElasticMatrix - diag(Mus);
Connectivities = sum(L>0);
NNodes = size(NodePositions,1);

[Edges,Lambdas,Mus] = DecodeElasticMatrix(ElasticMatrix);
NumberOfGraphs = 0;
for i=1:size(Edges,1)
    if(Connectivities(Edges(i,1))>1)&(Connectivities(Edges(i,2))>1)
        NumberOfGraphs=NumberOfGraphs+1;
    end
end

NodePositionArray = zeros(NNodes-1,size(NodePositions,2),NumberOfGraphs);
ElasticMatrices = zeros(NNodes-1,NNodes-1,NumberOfGraphs);
NodeIndicesArray = zeros(NNodes-1,NumberOfGraphs);

k=1;
for i=1:size(Edges,1)
    if(Connectivities(Edges(i,1))>1)&(Connectivities(Edges(i,2))>1)
        [em] = f_reattach_edges(ElasticMatrix,Edges(i,1),Edges(i,2));
        temp = NodePositions(Edges(i,2),:);
        [np,em,inds] = f_remove_node(NodePositions,em,Edges(i,2));
        np(Edges(i,1),:) = (np(Edges(i,1),:)+temp)/2;
        NodePositionArray(:,:,k) = np;
        ElasticMatrices(:,:,k) = em(:,:);
        NodeIndicesArray(:,k) = inds;
        k=k+1;
    end
end

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% some elementary graph transformations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [NodePositions2, ElasticMatrix2, NodeIndices] = f_remove_node(NodePositions,ElasticMatrix,NodeNumber)
% remove from the graph node number NodeNumber
    newinds = [1:NodeNumber-1,NodeNumber+1:size(NodePositions,1)];
    NodePositions2 = NodePositions(newinds,:);
    ElasticMatrix2 = ElasticMatrix(newinds,newinds);
    NodeIndices = newinds;
end

function [ElasticMatrix2] = f_reattach_edges(ElasticMatrix,NodeNumber1,NodeNumber2)
% reattaches all edges connected with NodeNumber2 to NodeNumber1
% and make a new star with an elasticity average of two merged stars
ElasticMatrix2(:,:) = ElasticMatrix;
mus = diag(ElasticMatrix);
lm = ElasticMatrix-diag(mus);
ElasticMatrix2(NodeNumber1,:) = max(lm(NodeNumber1,:),lm(NodeNumber2,:));
ElasticMatrix2(:,NodeNumber1) = max(lm(:,NodeNumber1),lm(:,NodeNumber2));
ElasticMatrix2(NodeNumber1,NodeNumber1) = (ElasticMatrix(NodeNumber1,NodeNumber1)+ElasticMatrix(NodeNumber2,NodeNumber2))/2;
end

function [NodePositions2, ElasticMatrix2, NodeIndices] =...
    f_add_nonconnected_node(NodePositions, ElasticMatrix, NewNodePosition)
% Add new node without connections
    NodePositions2(:,:) = NodePositions(:,:);
    NodePositions2(size(NodePositions,1)+1,:) = NewNodePosition(:); 
    ElasticMatrix2(:,:) = ElasticMatrix(:,:);
    ElasticMatrix2(size(NodePositions2,1),:) = zeros(size(NodePositions,1),1);
    ElasticMatrix2(:,size(NodePositions2,1)) = zeros(1,size(NodePositions2,1));
    NodeIndices = [1:size(NodePositions,1),0];
end

function [ElasticMatrix2] = f_removeedge(ElasticMatrix, Node1, Node2)
% remove edge connecting Node1 and Node 2
    ElasticMatrix2(:,:) = ElasticMatrix(:,:);
    ElasticMatrix2(Node1,Node2) = 0;
    ElasticMatrix2(Node2,Node1) = 0;
end

function [ElasticMatrix2] = f_add_edge(ElasticMatrix, Node1, Node2, lambda)
% connects Node1 and Node 2 by an edge with elasticity lambda
    ElasticMatrix2(:,:) = ElasticMatrix(:,:);
    ElasticMatrix2(Node1,Node2) = lambda;
    ElasticMatrix2(Node2,Node1) = lambda;
end

function [NodePositions2, ElasticMatrix2, NodeIndices] = f_get_star(NodePositions,ElasticMatrix,NodeCenter)
% extracts a star from the graph with the center in NodeCenter
NodeIndices = find(ElasticMatrix(NodeCenter,:)>0);
NodePositions2(:,:) = NodePositions(NodeIndices,:);
ElasticMatrix2(:,:) = ElasticMatrix(NodeIndices,NodeIndices);
end