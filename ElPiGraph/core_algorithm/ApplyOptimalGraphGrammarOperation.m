function [graphNew, partNew] = ...
    ApplyOptimalGraphGrammarOperation(data, graph, part, operationtypes)
% This functinon applies the most optimal graph grammar operation of
% operationtype for the embedment of an elastic graph described by
% ElasticMatrix.
%
%Inputs:
%   graph is structure of graph description and construction parameters.
%           Structure contains following fields:
%       nNodes is number of nodes;
%       NodePositions - is k-by-m matrix of positions of the graph nodes in
%           the same space as data. 
%       Lambda is predefined elasticity of edge
%       Mu is predefined elasticity of star
%       Lambdas is nNodes-by-nNodes symmetric matrix describing the
%           connectivity and the elastic properties of the graph.
%           Non-zero element Lambdas(i,j) is elasticity  modulo of edge
%           which connects nodes i and j. Main diagonal elements are zero.
%       Mus is vector of stars elasticities. Non-zero Mus(i) corresponds to
%       star with centre in node i.
%       TrimmingRadius is squared trimming radius, a parameter required for
%           robust principal graphs (see https://github.com/auranic/Elastic-principal-graphs/wiki/Robust-principal-graphs)
%       LocalSearch is boolean with value true corresponds to local search.
%       RadiusOfLocalSearch is graph distance to define neighbourhood for
%           local search.
%       MaxNumberOfIterations is maximal number of iteration for graph
%           embedding (usually two iterations is enough).
%       eps is minimal relative change in the node positions to be
%           considered the graph embedded (convergence criteria)
%       MaxMemorySize is the maximum memory size for the block of the
%           distance matrix when partition the data. This means that
%           maximal number of simultaneously processed data points is not
%           greater than MaxBlockSize/k, where k is number of nodes. 
%       MaxBlockSize is maximal number of simultaneously tested data
%           points: MaxBlockSize = MaxMemorySize / nNodes;
%   data is structure which contains all data related values. Structure
%           has following fields:
%       nPoints is number of data points.
%       dim is dimension (number of coordinates) of data
%       X is the nPoints-by-dim data matrix. Each row corresponds to one
%           data point.
%       Weights is nPoints-by-1 vector of data point weights (one weoght
%           for each point).
%       SquaredX is nPoints-by-1 vector of squared length of each data
%           vector.
%       XW = is the nPoints-by-dim matrix with data multiplied by data
%           point weights.
%   part is structure which contains data point associations with graph's
%           nodes. Structure contains following fields:
%       partition is nPoints-by-1 vector. partition(i) is number of node
%           which is associated with data point X(i,:). For robust version
%           of graph partition(i)==0 if point X(i,:) futher from all nodes
%           then TrimmingRadius.
%       dists is nPoints-by-1 vector. dists(i) is squared distance between
%           node with number partition(i) and data point X(i,:). For robust
%           version of graph dist(i) = TrimmingRadius^2 for all points with
%           partition(i)==0. 
%   operationtypes is cell array with names of graph grammar operations.
%
% Returns recalculated structures graph and part.

    % Check existance of partition
    if isempty(part.partition)
        part = PartitionDataInt(data, graph, part);
    end        

    % Provide vector with non-zero Mus if nNodes == 2
    if graph.nNodes == 2
        graph.Mus(:) = graph.Mu;
    end
    
    % Form new graphs to select the best
    NodePositionArrayAll = [];
    ElasticMatricesAll = [];
    ElasticVectorsAll = [];
    NodeIndicesArrayAll = [];
    for i=1:size(operationtypes,1)
        [NodePositionArray, ElasticMatrices,...
            ElsaticVectors, NodeIndicesArray] =... 
            GraphGrammarOperation(graph.NodePositions, graph.Lambdas,...
            graph.Mus, data.X, char(operationtypes(i)), part.partition);
        NodePositionArrayAll = cat(3, NodePositionArrayAll, NodePositionArray);
        ElasticMatricesAll   = cat(3, ElasticMatricesAll, ElasticMatrices);
        ElasticVectorsAll    = cat(2, ElasticVectorsAll, ElsaticVectors);
        NodeIndicesArrayAll  = cat(2, NodeIndicesArrayAll, NodeIndicesArray);
    end

    % special case - if the set of operations is empty then we optimize the
    % unmodified matrix
    
    if size(operationtypes,1)==0
        NodePositionArrayAll(:,:,1) = graph.NodePositions;
        ElasticMatricesAll(:, :, 1) = graph.Lambdas;
        ElasticVectorsAll(:, 1) = graph.Mus;
        NodeIndicesArrayAll = 1:size(graph.NodePositions,1);
    end
    
    
    % Currently found the best Energy
    minEnergy = realmax;
    %Just in case prepare array of indices for special case
    sd = graph.nNodes;
    sd1 = GetNeighbourhoodOnTheGraph(graph.Lambdas, sd, 1);
    sd1 = fast_setdiff1(sd1,sd);
    % Change some values of graph for new size
    graph.nNodes = size(NodePositionArrayAll, 1);
    graph.MaxBlockSize = floor(graph.MaxMemorySize / graph.nNodes);
    % Copy graphs and partitions
    graphNew = graph;
    partNew = part;
    
    
    % Test each possible continuation and select the best one
    for i=1:size(ElasticMatricesAll,3)
        % Form new graph to test.
        graph.NodePositions = NodePositionArrayAll(:, :, i);
        graph.Lambdas = ElasticMatricesAll(:, :, i);
        graph.Mus = ElasticVectorsAll(:, i);
        % Remove unnecessary non-zero mus if nNodes == 2
        if graph.nNodes == 3
            inds = sum(graph.Lambdas > 0) == 1;
            graph.Mus(inds) = 0;
        end
        if graph.LocalSearch
            NodeIndices = NodeIndicesArrayAll(:,i);
            sd = fast_setdiff1(1:graph.nNodes, NodeIndices);
            % Search the "initial" nodes.
            if isempty(sd)
                sd = sd1;
            end
            graph.NodesSubSet =...
                GetNeighbourhoodOnTheGraph(graph.Lambdas,...
                sd, graph.RadiusOfLocalSearch);
            [graph, part1, ElasticEnergy] =...
                PrimitiveElasticGraphEmbedment(data, graph, part);
        else 
            
            graph_temp = graph;
            if graph.PenalizedEnergy
                alpha = graph.BranchingControls(1);
                beta = graph.BranchingControls(2);
                Connectivity = sum(graph.Lambdas>0);
                graph_temp.Mus(Connectivity>2) = graph.Mus(Connectivity>2)/beta;
                stars = find(Connectivity>2);
                for i=1:length(stars)
                    leaves = find(graph_temp.Lambdas(:,stars(i))>0);
                    graph_temp.Lambdas(stars(i),leaves) = graph.Lambdas(stars(i),leaves)/beta;
                    graph_temp.Lambdas(leaves,stars(i)) = graph.Lambdas(leaves,stars(i))/beta;
                end
            end
            
            [graph_res, part1, ElasticEnergy] =...
                PrimitiveElasticGraphEmbedment(data, graph_temp, part);
            
            graph_res.Mus = graph.Mus;
            graph_res.Lambdas = graph.Lambdas;
            graph = graph_res;
        end
        % Penalizing energy (currently - control for branching, in the future
        % - it can be also penalizing empty nodes, number of graph components, etc.)
        if graph.PenalizedEnergy
            ElasticEnergy = ComputePenalizedElasticEnergy(data,graph,part1);
        end
        % Remember the best
        if ElasticEnergy<minEnergy
            minEnergy = ElasticEnergy;
            graphNew = graph;
            partNew = part1;
        end
    end
end

function Z = fast_setdiff1(X,Y)
    if ~isempty(X) && ~isempty(Y)
        X = X + 1;
        Y = Y + 1;
        check = false(1, max(max(X), max(Y)));
        check(X) = true;
        check(Y) = false;
        Z = X(check(X));
        Z = Z - 1;
    else
        Z = X;
    end
end

