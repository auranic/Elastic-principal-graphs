function [NodePositions2, ElasticMatrix2, partition, dists] = ...
    ApplyOptimalGraphGrammarOpeation(X, NodePositions, ElasticMatrix,...
    operationtypes, varargin)
% This functinon applies the most optimal graph grammar operation of
% operationtype for the embedment of an elastic graph described by
% ElasticMatrix.
%
%Inputs:
%   X - is the n-by-m data matrix. Each row corresponds to one data point.
%   NodePositions - is k-by-m matrix of positions of the original graph
%       nodes in the same space as X.
%   ElasticMatrix - k-by-k symmetric matrix describing the connectivity and
%       the elastic properties of the graph. Star elasticities (mu
%       coefficients) are presented on the main diagonal (non-zero entries
%       only for star centers), and the edge elasticity moduli are
%       presented out of diagonal.
%   operationtypes is cell array with names of graph grammar operations.
%   varargin contains Name, Value pairs. Names can be: 
%   'MaxBlockSize' with integer number which is maximum size of the block
%       of the distance matrix when partition the data. This means that
%       maximal size of distance matrix is MaxBlockSize-by-k where k is
%       number of nodes.
%   'verbose' with 1/0 is to display/hide the energy values at each
%       iteration and in the end of the process.
%   'TrimmingRadius' is trimming radius, a parameter required for robust
%       principal graphs (see https://github.com/auranic/Elastic-principal-graphs/wiki/Robust-principal-graphs)
%   'LocalSearch' specifies local test for optimizing only a subset of the
%       nodes. Integer number defines radius of neghbourhood to use.
%   The following Name,Value pairs can be used in called functions
%   'MaxNumberOfIterations' with integer number which is maximum number of
%       iterations for EM algorithm. 
%   'eps' with real number which is minimal relative change in the node
%       positions to be considered the graph embedded (convergence criteria) 
%   'PointWeights' with n-by-1 vector of data point weights
%
% Outputs
%   NodePositions2 is positions of empbedded nodes.
%   ElasticMatrix2 is elastic matrix of the best selected new graph.
%   partition is n-by-1 vector. partition(i) is number of node which is
%       associated with data point X(i,:).
%   dists is array of squared distances form each data point to nerest
%       node.

    % Parse input arguments
    LocalSearch = 0;
    RadiusOfLocalSearch = 0;
    MaxBlockSize = 100000;
    TrimmingRadius = Inf;

    for i=1:2:length(varargin)
        if strcmpi(varargin{i},'LocalSearch')
            LocalSearch = 1;
            RadiusOfLocalSearch = varargin{i+1};
        end
        if strcmpi(varargin{i},'MaxBlockSize')
            MaxBlockSize = varargin{i+1};
        end
        if strcmpi(varargin{i},'TrimmingRadius')
            TrimmingRadius = varargin{i+1};
        end
    end
    TrimmingRadius = TrimmingRadius .^ 2;

    % We compute these things here in order not to recompute for each graph
    % embedment 
    SquaredX = sum(X.^2, 2);
    partition = PartitionData(X, NodePositions, MaxBlockSize,...
        SquaredX, TrimmingRadius);

    % Form new graphs to select the best
    NodePositionArrayAll = [];
    ElasticMatricesAll = [];
    NodeIndicesArrayAll = [];
    for i=1:size(operationtypes,1)
        [NodePositionArray, ElasticMatrices, NodeIndicesArray] =... 
            GraphGrammarOperation(NodePositions, ElasticMatrix,...
                X, char(operationtypes(i)), partition);
        NodePositionArrayAll = cat(3, NodePositionArrayAll, NodePositionArray);
        ElasticMatricesAll = cat(3, ElasticMatricesAll, ElasticMatrices);
        NodeIndicesArrayAll = cat(2, NodeIndicesArrayAll, NodeIndicesArray);
    end

    minEnergy = realmax;

    % Form data for local search
    if LocalSearch
        LocalInfo = struct();
        LocalInfo.Partition = [];
        [LocalInfo.Partition] = partition;
    end

    % Test each possible continuation and select the best one
    for i=1:size(ElasticMatricesAll,3)
        %Control of ElasticMatrix
        em = ElasticMatricesAll(:,:,i);
        mu = diag(em);
        em1 = em - diag(mu);
        inds = sum(em1 > 0) == 1;
        mu(inds) = 0;
        em = em1 + diag(mu);
        if LocalSearch
            NodeIndices = NodeIndicesArrayAll(:,i);
            sd = fast_setdiff1(1:size(NodePositionArrayAll(:, :, i), 1),...
                NodeIndices);
            if size(sd, 2) == 0
                sd = size(NodePositions, 1);
                sd1 = GetNeighbourhoodOnTheGraph(ElasticMatrix, sd, 1);
                sd = fast_setdiff1(sd1,sd);
            end
            LocalInfo.Nodes =...
                GetNeighbourhoodOnTheGraph(ElasticMatricesAll(:, :, i),...
                sd, RadiusOfLocalSearch);
            [np, ElasticEnergy, part, dist] =...
                PrimitiveElasticGraphEmbedment(X,...
                NodePositionArrayAll(:, :, i), em,...
                'verbose', 0, 'SquaredX', SquaredX, 'Local', LocalInfo);
        else
            [np, ElasticEnergy, part, dist] =...
                PrimitiveElasticGraphEmbedment(X,...
                NodePositionArrayAll(:,:,i), em,...
                'verbose', 0, 'SquaredX', SquaredX);
        end
        % Remember the best
        if ElasticEnergy<minEnergy
            minEnergy = ElasticEnergy;
            NodePositions2 = np;
            partition = part;
            dists = dist;
            ElasticMatrix2 = em;
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

