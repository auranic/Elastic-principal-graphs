function [EmbeddedNodePositions, ElasticEnergy, partition, MSE,EP,RP] =...
    PrimitiveElasticGraphEmbedment(X, NodePositions, ElasticMatrix, varargin)
% This is the core function fitting a primitive elastic graph to the data
%Inputs
%   X - is the n-by-m data matrix.
%   NodePositions - is k-by-m matrix of positions of the graph nodes in the
%       same space as X.
%   ElasticMatrix - k-by-k a matrix describing the connectivity of the
%       graph. Star elasticities (mu coefficients) are presented along the
%       main diagonal (non-zero entries only for star centers), and the
%       weighted adjacency matrix at non-diagonal elements.
%   varargin contains Name, Value pairs. Names can be: 
%   'MaxNumberOfIterations' with integer number which is maximum number of
%       iterations for EM algorithm. 
%   'eps' with real number which is minimal relative change in the node
%       positions to be considered the graph embedded (convergence criteria) 
%   'PointWeights' with n-by-1 vector of data point weights
%   'MaxBlockSize' with integer number which is maximum size of the block
%       of the distance matrix when partition the data. This means that
%       maximal size of distance matrix is MaxBlockSize-by-k where k is
%       number of nodes.
%   'verbose' with 1/0 is to display/hide the energy values at each
%       iteration and in the end of the process.
%   'TrimmingRadius' is trimming radius (link ???)
%

    MaxNumberOfIterations = 10;
    eps = 0.01;
    verbose = 0;
    TrimmingRadius = Inf;

    N = size(X,1);
    MaxBlockSize = 100000;

    PointWeights = ones(N,1);

    for i=1:2:length(varargin)
        if strcmpi(varargin{i},'MaxNumberOfIterations')
            MaxNumberOfIterations = varargin{i+1};
            if isdeployed
                MaxNumberOfIterations = str2double(MaxNumberOfIterations);
            end
        elseif strcmpi(varargin{i},'eps')
            eps = varargin{i+1};
            if isdeployed
                eps = str2double(eps);
            end
        elseif strcmpi(varargin{i},'MaxBlockSize')
            MaxBlockSize = varargin{i+1};
            if isdeployed
                MaxBlockSize = str2double(MaxBlockSize);
            end
        elseif strcmpi(varargin{i},'verbose')
            verbose = varargin{i+1};
            if isdeployed
                verbose = str2double(verbose);
            end
        elseif strcmpi(varargin{i},'TrimmingRadius')
            TrimmingRadius = varargin{i+1};
            if isdeployed
                TrimmingRadius = str2double(TrimmingRadius);
            end
        elseif strcmpi(varargin{i},'PointWeights')
            PointWeights = varargin{i+1};
        end
    end
    
    %Auhiliary computations
    [SpringLaplacianMatrix] = ComputeSpringLaplacianMatrix(ElasticMatrix);
    SquaredX = sum(X.^2, 2);

    % Main iterative EM cycle: partition, fit given the partition, repeat
    for i = 1:MaxNumberOfIterations
        if verbose
            [partition, dists] = PartitionData(X, NodePositions, MaxBlockSize, SquaredX, TrimmingRadius);
            [ElasticEnergy, MSE, EP, RP] =...
                ComputePrimitiveGraphElasticEnergy(X, NodePositions, ElasticMatrix, partition, dists, 0);
        else
            [partition] = PartitionData(X, NodePositions, MaxBlockSize, SquaredX, TrimmingRadius);
            ElasticEnergy = 0;
        end

        [NewNodePositions] =...
            FitGraph2DataGivenPartition(X, PointWeights, NodePositions, SpringLaplacianMatrix, partition);

        diff = ComputeRelativeChangeOfNodePositions(NodePositions, NewNodePositions);

        if verbose
            display(sprintf('Iteration %i, diff=%3.3f, E=%4.4f, MSE=%4.4f, EP=%4.4f, RP=%4.4f,',...
                i, diff, ElasticEnergy, MSE, EP, RP));
        end

        if(diff<eps) break; end;

        NodePositions = NewNodePositions;
    end

    [partition, dists] = PartitionData(X, NodePositions, MaxBlockSize, SquaredX, TrimmingRadius);
    [ElasticEnergy, MSE, EP, RP] =...
        ComputePrimitiveGraphElasticEnergy(X, NodePositions, ElasticMatrix, partition, dists, 0);
    if verbose
        display(sprintf('E=%4.4f, MSE=%4.4f, EP=%4.4f, RP=%4.4f,', ElasticEnergy, MSE, EP, RP));
    end
    EmbeddedNodePositions = NodePositions;
end


%%%%%%%%%%%%%%%%%%%%%%%
%% Transforms the ElasticMatrix into the SpringLaplacianMatrix ready to be used in the SLAU solving
%%%%%%%%%%%%%%%%%%%%%%%
function [SpringLaplacianMatrix] = ComputeSpringLaplacianMatrix(ElasticMatrix)

NumberOfNodes = size(ElasticMatrix,1);

% first, make the vector of mu coefficients
Mu = diag(ElasticMatrix);
% create the weighted adjacency matrix
Lambda = ElasticMatrix - diag(Mu);
% Diagonal matrix of weighted connectivities
LambdaSums = sum(Lambda);
DL = diag(LambdaSums);
% E matrix (contribution from edges) is simply weighted Laplacian
E = DL-Lambda;

% S matrix (contribution from stars) is composed of Laplacian for positive strings (star edges) with
% elasticities mu/k, where k is the order of the star, and Laplacian for
% negative strings with elasticities -mu/k^2. Negative springs connect all
% star leaves in a clique.

StarCenterIndices = find(Mu>0);
S = zeros(NumberOfNodes,NumberOfNodes);

for i=1:size(StarCenterIndices,1)
    Spart = zeros(NumberOfNodes,NumberOfNodes);
        % leaf indices
        leafs = find(Lambda(:,StarCenterIndices(i))>0);
        % order of the star
        K = size(leafs,1);
        
        Spart(StarCenterIndices(i),StarCenterIndices(i)) = Mu(StarCenterIndices(i));
        Spart(StarCenterIndices(i),leafs) = -Mu(StarCenterIndices(i))/K;
        Spart(leafs,StarCenterIndices(i)) = -Mu(StarCenterIndices(i))/K;
        Spart(leafs,leafs) = Mu(StarCenterIndices(i))/K^2;
    S = S+Spart;
end

SpringLaplacianMatrix = E+S;
end

function [NodeClusterCenters, NodeClusterRelativeSize] =...
    ComputeWeightedAverage(X, partition, PointWeights, NumberOfNodes)
%ComputeWeightedAverage calculate NodeClausterCentres as weighted averages
%of points from matrix X.
%
%Inputs
%   X is n-by-m matrix of data points where each row corresponds to one
%       observation.
%   partition is n-by-1 (column) vector of node numbers. This vector
%       associate data points with Nodes.
%   PointWeights is n-by-m (column) vector of point weights.
%   NumberOfNodes is number of nodes to calculate means.
%
%Important! if there is no point associated with node then coordinates of
%this node centroid are zero.
%
     NodeClusterCenters = zeros(NumberOfNodes,size(X,2));
     NodeClusterRelativeSize = zeros(NumberOfNodes,1);
     TotalWeight = sum(PointWeights);
     for i=1:NumberOfNodes
         inds = find(partition==i);
         if(size(inds,1)>0)
             tmp = sum(PointWeights(inds));
             NodeClusterCenters(i,:) = sum(bsxfun(@times,X(inds,:),PointWeights(inds)),1)/tmp;
             NodeClusterRelativeSize(i) = tmp/TotalWeight;
         end
     end
end



%%%%%%%%%%%%%%%%%%%%%%%
%% Solves the SLAU to find new node positions
%%%%%%%%%%%%%%%%%%%%%%%
function [NewNodePositions] =...
    FitGraph2DataGivenPartition(X, PointWeights, NodePositions,...
    SpringLaplacianMatrix, partition)

    NumberOfNodes = size(SpringLaplacianMatrix,1);
    [NodeClusterCenters, NodeClusterRelativeSize] =...
        ComputeWeightedAverage(X, partition, PointWeights,NumberOfNodes);
    SLAUMatrix = diag(NodeClusterRelativeSize)+SpringLaplacianMatrix;
    NewNodePositions = SLAUMatrix\bsxfun(@times, NodeClusterRelativeSize, NodeClusterCenters);
end

%%%%%%%%%%%%%%%%%%%%%%%
%% Estimates the relative difference between two node configurations
%%%%%%%%%%%%%%%%%%%%%%%
function [diff] = ComputeRelativeChangeOfNodePositions(NodePositions,NewNodePositions)
    %diff = norm(NodePositions-NewNodePositions)/norm(NewNodePositions);
    %diff = immse(NodePositions,NewNodePositions)/norm(NewNodePositions);
    diff = sum((NodePositions-NewNodePositions).^2,2)./sum(NewNodePositions.^2,2);
    diff = max(diff);
    %diff = mean(diff);
end