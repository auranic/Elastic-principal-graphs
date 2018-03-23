function [graph, part, ElasticEnergy, MSE, EP, RP] =...
    PrimitiveElasticGraphEmbedment(data, graph, part)
% This is the core function for fitting a primitive elastic graph to the data
% Inputs
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
%       NodesSubSet is subset of nodes to use in local search.
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
%
% Outputs
%   graph is structure with embegged graph
%   part is structure with the results of the last partition
%   ElasticEnergy is total elastic energy 
%   MSE is mean square error of data approximation.
%   EP is edge potential 
%   RP is harmonicity potential 
%

    %Auxiliary computations    
    SpringLaplacianMatrix = ComputeSpringLaplacianMatrix(graph);

    loc = graph.LocalSearch;
    if graph.nNodes < 6
        graph.LocalSearch = false;
    end
    
    if ~graph.LocalSearch
        % Main iterative EM cycle: partition, fit given the partition, repeat
        for i = 1:graph.MaxNumberOfIterations
            part = PartitionDataInt(data, graph, part);
            % Fit new graph
            NewNodePositions =...
                FitGraph2DataGivenPartition(data, part, SpringLaplacianMatrix);
            % Calculate size of change
            diff = ComputeRelativeChangeOfNodePositions(graph.NodePositions,...
                NewNodePositions);
            % Save new result
            graph.NodePositions = NewNodePositions;
            % Stop is cahnge is small enough
            if diff < graph.eps
                break; 
            end;
        end
    else
        % Local data is the data that does not belong to the fixed nodes
        % LocalInds is logical index with value true for data points which
        % does not correspond to complement of selected set of nodes
        % NodeSubSet.
        % NodeSubSet is array with numbers of nodes to position optimisation.

        % Form complement set
        FixedSubSet = 1:graph.nNodes;
        FixedSubSet(graph.NodesSubSet) = [];
        % Define set of points to use
        LocalInds = ~ismember(part.partition, FixedSubSet);
        % Extract corresponding part of arrays
        dataLocal.X = data.X(LocalInds,:);
        dataLocal.nPoints = sum(LocalInds);
        dataLocal.dim = data.dim;
        dataLocal.Weights = data.Weights(LocalInds);
        dataLocal.SquaredX = data.SquaredX(LocalInds);
        dataLocal.XW = data.XW(LocalInds,:);
        partLocal.partition = part.partition(LocalInds);
        partLocal.dists = part.dists(LocalInds);

        for i=1:graph.MaxNumberOfIterations
            partLocal = PartitionDataInt(dataLocal, graph,...
                partLocal, graph.NodesSubSet);
            % Fit new graph
            NewNodePositions =...
                FitGraph2DataGivenPartitionLocal(dataLocal, graph,...
                partLocal, SpringLaplacianMatrix);
            % Calculate size of change
            diff = ComputeRelativeChangeOfNodePositions(...
                graph.NodePositions,...
                NewNodePositions);
            % Save new result
            graph.NodePositions = NewNodePositions;
            % Stop is cahnge is small enough
            if diff < graph.eps
                break; 
            end;
        end
    end
    % Restore saved value
    graph.LocalSearch = loc;
    
    % Recalculate partition for the final state
    part = PartitionDataInt(data, graph, part);
    % Compute energies
    [ElasticEnergy, MSE, EP, RP] =...
        ComputePrimitiveGraphElasticEnergy(data, graph, part);
end

function SpringLaplacianMatrix =...
    ComputeSpringLaplacianMatrix(graph)
%%%%%%%%%%%%%%%%%%%%%%%
%% Transforms the ElasticMatrix into the SpringLaplacianMatrix ready 
%% to be used in the SLAU solving 
%%%%%%%%%%%%%%%%%%%%%%%

    % Diagonal matrix of edge elasticities
    LambdaSums = sum(graph.Lambdas);
    % E matrix (contribution from edges) is simply weighted Laplacian
    E = diag(LambdaSums) - graph.Lambdas;

    % matrix S (contribution from stars) is composed of Laplacian for
    % positive strings (star edges) with elasticities mu/k, where k is the
    % order of the star, and Laplacian for negative strings with
    % elasticities -mu/k^2. Negative springs connect all star leafs in a
    % clique. 

    StarCenterIndices = find(graph.Mus>0);
    
    S = zeros(graph.nNodes, graph.nNodes);
    
    for i=1:size(StarCenterIndices,1)
        Spart = zeros(graph.nNodes,graph.nNodes);
        % leaf indices
        leafs = graph.Lambdas(:,StarCenterIndices(i))>0;
        % order of the star
        K = sum(leafs);
        
        Spart(StarCenterIndices(i),StarCenterIndices(i)) =...
            graph.Mus(StarCenterIndices(i));
        Spart(StarCenterIndices(i),leafs) = -graph.Mus(StarCenterIndices(i)) / K;
        Spart(leafs,StarCenterIndices(i)) = -graph.Mus(StarCenterIndices(i)) / K;
        Spart(leafs,leafs) = graph.Mus(StarCenterIndices(i))/K^2;
        S = S + Spart;
    end
    SpringLaplacianMatrix = E + S;
end

function [NodeClusterCenters, NodeClusterRelativeSize] =...
    ComputeWeightedAverage(data, part, NumberOfNodes)
%ComputeWeightedAverage calculate NodeClusterCentres as weighted averages
%of points from matrix X.
%
%Important! if there is no point associated with node then coordinates of
%this node centroid are zero.
%
    % Auxiliary calculations
    partp1 = part.partition + 1;
    % Calculate total weights
    TotalWeight = sum(data.Weights);
    % Calculate weights for Relative size
    tmp = accumarray(partp1, data.Weights, [NumberOfNodes + 1, 1]);
    NodeClusterRelativeSize = tmp(2:end) / TotalWeight;
    % To prevent appearance of NaN
    tmp(tmp == 0) = 1;
    NodeClusterCenters = zeros(NumberOfNodes + 1, data.dim);
    for k=1:data.dim
        NodeClusterCenters(:, k) =...
            accumarray(partp1, data.XW(:, k), [NumberOfNodes+1, 1]) ./ tmp;
    end
    NodeClusterCenters = NodeClusterCenters(2:end,:);
end

function NewNodePositions =...
    FitGraph2DataGivenPartition(data, part, SpringLaplacianMatrix)
%%%%%%%%%%%%%%%%%%%%%%%
%% Solves the SLAU to find new node positions
%%%%%%%%%%%%%%%%%%%%%%%

    NumberOfNodes = size(SpringLaplacianMatrix, 1);
    [NodeClusterCenters, NodeClusterRelativeSize] =...
        ComputeWeightedAverage(data, part, NumberOfNodes);
    SLAUMatrix = diag(NodeClusterRelativeSize) + SpringLaplacianMatrix;
    NewNodePositions = SLAUMatrix...
        \bsxfun(@times, NodeClusterRelativeSize, NodeClusterCenters);
end

function [NewNodePositions] =...
    FitGraph2DataGivenPartitionLocal(data, graph, part, SpringLaplacianMatrix)

% XLocal, PointWeightsLocal,...
%     NodePositions, SpringLaplacianMatrix, partitionLocal, NodeSubSet)
%%%%%%%%%%%%%%%%%%%%%%%
%% Solves the SLAU to find new node positions, local version
%%%%%%%%%%%%%%%%%%%%%%%
    % Define sizes
%     NumberOfNodes = graph.nNodes;
    SizeSubSet = length(graph.NodesSubSet);
    [NodeClusterCentersLocal, NodeClusterRelativeSizeLocal] =...
        ComputeWeightedAverage(data, part, SizeSubSet);

    rs = zeros(graph.nNodes, 1);
    rs(graph.NodesSubSet, :) = NodeClusterRelativeSizeLocal(:,:);
    SLAUMatrix = diag(rs)+SpringLaplacianMatrix;
    
    SLAUMatrixLocal = SLAUMatrix(graph.NodesSubSet,graph.NodesSubSet);
    
    ComplementNodeSet = 1:graph.nNodes;
    ComplementNodeSet(graph.NodesSubSet) = [];
    
    ComplementSLAUMatrix = SLAUMatrix(graph.NodesSubSet,ComplementNodeSet);
    
    RightHandSide = bsxfun(@times, NodeClusterRelativeSizeLocal,...
        NodeClusterCentersLocal);
    
    rhs1 = graph.NodePositions(ComplementNodeSet, :);
    
    RightHandSide = RightHandSide - ComplementSLAUMatrix * rhs1;
    
    NewNodePositionsLocal = SLAUMatrixLocal \ RightHandSide;
    
    NewNodePositions(ComplementNodeSet, :) =...
        graph.NodePositions(ComplementNodeSet, :);
    NewNodePositions(graph.NodesSubSet, :) = NewNodePositionsLocal;
end

function diff =...
    ComputeRelativeChangeOfNodePositions(NodePositions, NewNodePositions)
    diff = max(sum((NodePositions - NewNodePositions) .^ 2, 2)...
        ./ sum(NewNodePositions .^ 2, 2));
end