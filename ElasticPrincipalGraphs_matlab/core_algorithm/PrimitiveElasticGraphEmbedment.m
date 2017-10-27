function [EmbeddedNodePositions, ElasticEnergy, partition, MSE,EP,RP]...
    = PrimitiveElasticGraphEmbedment(X, NodePositions, ElasticMatrix, varargin)
% This is the core function for fitting a primitive elastic graph to the data
% Inputs
%   X - is the n-by-m data matrix. Each row corresponds to one data point.
%   NodePositions - is k-by-m matrix of positions of the graph nodes in the
%       same space as X.
%   ElasticMatrix - k-by-k symmetric matrix describing the connectivity and
%       the elastic properties of the graph. Star elasticities (mu
%       coefficients) are presented on the main diagonal (non-zero entries
%       only for star centers), and the edge elasticity moduli are
%       presented out of diagonal.  
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
%   'TrimmingRadius' is trimming radius, a parameter required for robust
%       principal graphs (see https://github.com/auranic/Elastic-principal-graphs/wiki/Robust-principal-graphs)
%   'Local' specifies local test for optimizing only a subset of the nodes
%       In this case, the argument contains a structure with the following
%       fields 'Nodes' - which nodes are optimized, 'Partition' -
%       pre-computed partitioning of X
%   'SquaredX' with n-by-1 vector of squared length of data vectors.

    % Set profile = 1 for debug purposes
    profile = 0;

   if profile == 1
        global numberOfLocalPoints;
        global numberOfFits;
        global numberOfGraphNodes;
        global FractionOfGraphNodes;
        global TimeForFitting;
        global CountNumberOfIterations;
        global numberOfFullFits;
   end


    MaxNumberOfIterations = 10;
    eps = 0.01;
    verbose = 0;
    TrimmingRadius = Inf;

    N = size(X,1); %Number of data points.
    MaxBlockSize = 100000;

    PointWeights = ones(N,1);

    localVersion = 0;
    NodeSubSet = [];
    partition = [];
    SquaredX = [];

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
        elseif strcmpi(varargin{i},'Local')
            LocalInfo = varargin{i+1};
            localVersion = 1;
            NodeSubSet = LocalInfo.Nodes;
            partition = LocalInfo.Partition;
        %% version 1.1    
        elseif strcmpi(varargin{i},'SquaredX') 
            SquaredX = varargin{i+1};
        end
    end

    %Auxiliary computations    
    SpringLaplacianMatrix = ComputeSpringLaplacianMatrix(ElasticMatrix);
    if size(SquaredX ,1)==0
        SquaredX = sum(X .^ 2, 2);
    end

    if ~localVersion
        
        % Main iterative EM cycle: partition, fit given the partition, repeat
        for i=1:MaxNumberOfIterations
            
            if profile==1
                tic;
            end
            
            if verbose
                [partition, dists] = PartitionData(X, NodePositions, MaxBlockSize, SquaredX, TrimmingRadius);
                [ElasticEnergy, MSE, EP, RP] =...
                    ComputePrimitiveGraphElasticEnergy(X, NodePositions, ElasticMatrix, partition, dists, 0);
            else
                [partition] = PartitionData(X, NodePositions, MaxBlockSize, SquaredX, TrimmingRadius);
                ElasticEnergy = 0;
            end
            
            if profile == 1
                numberOfLocalPoints(end + 1) = N; %#ok<AGROW>
                numberOfFits = numberOfFits + 1;
                numberOfFullFits = numberOfFullFits + 1;
                numberOfGraphNodes = numberOfGraphNodes + size(NodePositions, 1);
                FractionOfGraphNodes(end + 1) = 1; %#ok<AGROW>
                TimeForFitting(end + 1) = toc; %#ok<AGROW>
            end
            
            [NewNodePositions] =...
                FitGraph2DataGivenPartition(X, PointWeights,...
                SpringLaplacianMatrix, partition);
            
            diff = ComputeRelativeChangeOfNodePositions(NodePositions,NewNodePositions);
            
            if verbose
                display(sprintf(['Iteration %i, diff=%3.3f, E=%4.4f,',...
                    ' MSE=%4.4f, EP=%4.4f, RP=%4.4f,'],...
                    i, diff, ElasticEnergy, MSE, EP, RP));
            end
            
            if(diff<eps) 
                break; 
            end;
            
            NodePositions = NewNodePositions;
        end
    else
        % now, local version when only a NodeSubSet is optimized
        %% version 1.1
        % Local data is the data that does not belong to the fixed nodes
        
        
%         %Debugging
%         fprintf('Number of nodes %d number of nodes to modify %d\n',...
%             size(NodePositions, 1), length(NodeSubSet));
%         
%         % we have to define it through negation
%         % because it allows using pre-computed partitioning
%         % otherwise we would need to know the partitioning for the optimized
%         % nodes also (and recompute the full partitioning for each tested graph)
%         inds = [];
%         AllPointIndices = 1:N; 
%         AllNodeIndices = 1:size(NodePositions, 1);
%         FixedNodeIndices = fast_setdiff(AllNodeIndices, NodeSubSet);
%         SizeFixedNodes = size(FixedNodeIndices, 2);
%         for i=1:SizeFixedNodes
%             inds = [inds;find(partition==FixedNodeIndices(i))];
%         end
%         LocalInds = fast_setdiff(AllPointIndices,inds);
%         XLocal = X(LocalInds,:);
%         PointWeightsLocal = PointWeights(LocalInds,:);
%         SquaredXLocal = SquaredX(LocalInds); %sum(XLocal.^2, 2);

        % LocalInds is logical index with true for data points which does
        % not correspond to complement of selected set of nodes NodeSubSet
        % Assumption: NodeSubSet is array with numbers of nodes to position
        % optimisation.
        % Form complement set
        FixedSubSet = 1 : size(NodePositions, 1);
        FixedSubSet(NodeSubSet) = [];
        % Define set of points to use
        LocalInds = ~ismember(partition, FixedSubSet);
        % Extract corresponding part of arrays
        XLocal = X(LocalInds,:);
        PointWeightsLocal = PointWeights(LocalInds,:);
        SquaredXLocal = SquaredX(LocalInds); %sum(XLocal.^2, 2);
        

        
        for i=1:MaxNumberOfIterations
            if profile == 1
                tic;
            end
            if verbose
                [partitionLocal, distsLocal] = PartitionData(XLocal,...
                    NodePositions(NodeSubSet,:), MaxBlockSize,...
                    SquaredXLocal, TrimmingRadius);
                [ElasticEnergy, MSE, EP, RP] =...
                    ComputePrimitiveGraphElasticEnergy(XLocal,...
                    NodePositions, ElasticMatrix,...
                    partitionLocal, distsLocal, 0);
            else
                [partitionLocal] = PartitionData(XLocal,...
                    NodePositions(NodeSubSet,:),...
                    MaxBlockSize, SquaredXLocal, TrimmingRadius);
                ElasticEnergy = 0;
            end
            
            [NewNodePositions] =...
                FitGraph2DataGivenPartitionLocal(XLocal,...
                PointWeightsLocal, NodePositions, SpringLaplacianMatrix,...
                partitionLocal, NodeSubSet);
            
            if profile == 1
                numberOfLocalPoints(end + 1) = size(XLocal,1); %#ok<AGROW>
                numberOfFits = numberOfFits + 1;
                numberOfGraphNodes = numberOfGraphNodes + size(NodeSubSet, 1);
                FractionOfGraphNodes(end + 1) = size(NodeSubSet, 1)...
                    / size(NodePositions, 1); %#ok<AGROW>
                TimeForFitting(end + 1) = toc; %#ok<AGROW>
            end
            
            diff = ComputeRelativeChangeOfNodePositions(NodePositions(NodeSubSet,:),NewNodePositions(NodeSubSet,:));
            
            if verbose
                display(sprintf(['Iteration %i, diff=%3.3f, E=%4.4f,',...
                    ' MSE=%4.4f, EP=%4.4f, RP=%4.4f,'],...
                    i, diff, ElasticEnergy, MSE, EP, RP));
            end
            
            if(diff<eps) 
                break; 
            end;
            
            NodePositions = NewNodePositions;
        end
        
    end

    if profile == 1
        if MaxNumberOfIterations>0
            CountNumberOfIterations(end + 1) = i;
        end
    end
     
    [partition, dists] = PartitionData(X, NodePositions, MaxBlockSize, SquaredX, TrimmingRadius);
    [ElasticEnergy, MSE, EP, RP] =...
        ComputePrimitiveGraphElasticEnergy(X, NodePositions, ElasticMatrix, partition, dists, 0);
    if verbose
        display(sprintf('E=%4.4f, MSE=%4.4f, EP=%4.4f, RP=%4.4f,', ElasticEnergy, MSE, EP, RP));
    end
    EmbeddedNodePositions = NodePositions;
end


function [SpringLaplacianMatrix] = ComputeSpringLaplacianMatrix(ElasticMatrix)
%%%%%%%%%%%%%%%%%%%%%%%
%% Transforms the ElasticMatrix into the SpringLaplacianMatrix ready to be used in the SLAU solving
%%%%%%%%%%%%%%%%%%%%%%%

NumberOfNodes = size(ElasticMatrix,1);

% first, make the vector of mu coefficients
Mu = diag(ElasticMatrix);
% create the matrix with edge elasticity moduli at non-diagonal elements
Lambda = ElasticMatrix - diag(Mu);
% Diagonal matrix of edge elasticities
LambdaSums = sum(Lambda);
DL = diag(LambdaSums);
% E matrix (contribution from edges) is simply weighted Laplacian
E = DL-Lambda;

% S matrix (contribution from stars) is composed of Laplacian for positive strings (star edges) with
% elasticities mu/k, where k is the order of the star, and Laplacian for
% negative strings with elasticities -mu/k^2. Negative springs connect all
% star leafs in a clique.

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

function NewNodePositions =...
    FitGraph2DataGivenPartition(X, PointWeights,...
    SpringLaplacianMatrix, partition)
%%%%%%%%%%%%%%%%%%%%%%%
%% Solves the SLAU to find new node positions
%%%%%%%%%%%%%%%%%%%%%%%

    NumberOfNodes = size(SpringLaplacianMatrix, 1);
    [NodeClusterCenters, NodeClusterRelativeSize] =...
        ComputeWeightedAverage(X, partition, PointWeights, NumberOfNodes);
    SLAUMatrix = diag(NodeClusterRelativeSize)+SpringLaplacianMatrix;
    NewNodePositions = SLAUMatrix...
        \bsxfun(@times, NodeClusterRelativeSize, NodeClusterCenters);
end



function [NewNodePositions] =...
    FitGraph2DataGivenPartitionLocal(XLocal, PointWeightsLocal,...
    NodePositions, SpringLaplacianMatrix, partitionLocal, NodeSubSet)
%%%%%%%%%%%%%%%%%%%%%%%
%% Solves the SLAU to find new node positions, local version
%%%%%%%%%%%%%%%%%%%%%%%

    NumberOfNodes = size(SpringLaplacianMatrix,1);
    SizeSubSet = max(size(NodeSubSet,1),size(NodeSubSet,2));
    [NodeClusterCentersLocal, NodeClusterRelativeSizeLocal] =...
        ComputeWeightedAverage(XLocal, partitionLocal, PointWeightsLocal, SizeSubSet);

    rs = zeros(size(SpringLaplacianMatrix,1),1);
    rs(NodeSubSet,:) = NodeClusterRelativeSizeLocal(:,:);
    SLAUMatrix = diag(rs)+SpringLaplacianMatrix;
    
    SLAUMatrixLocal = SLAUMatrix(NodeSubSet,NodeSubSet);
    
    ComplementNodeSet = fast_setdiff(1:NumberOfNodes,NodeSubSet);
    
    ComplementSLAUMatrix = SLAUMatrix(NodeSubSet,ComplementNodeSet);
    
    RightHandSide = bsxfun(@times, NodeClusterRelativeSizeLocal, NodeClusterCentersLocal);
    
    rhs1 = NodePositions(ComplementNodeSet,:);
    
    RightHandSide = RightHandSide - ComplementSLAUMatrix*rhs1;
    
    NewNodePositionsLocal = SLAUMatrixLocal\RightHandSide;
    
    NewNodePositions(ComplementNodeSet,:) = NodePositions(ComplementNodeSet,:);
    NewNodePositions(NodeSubSet,:) = NewNodePositionsLocal(:,:);
end

function [diff] = ComputeRelativeChangeOfNodePositions(NodePositions,NewNodePositions)
%%%%%%%%%%%%%%%%%%%%%%%
%% Estimates the relative difference between two node configurations
%%%%%%%%%%%%%%%%%%%%%%%

    %diff = norm(NodePositions-NewNodePositions)/norm(NewNodePositions);
    %diff = immse(NodePositions,NewNodePositions)/norm(NewNodePositions);
    diff = sum((NodePositions-NewNodePositions).^2,2)./sum(NewNodePositions.^2,2);
    diff = max(diff);
    %diff = mean(diff);
end

%% faster alternative for setdiff
function [Z] = fast_setdiff(X,Y)
  check = false(1, max(max(X), max(Y)));
  check(X) = true;
  check(Y) = false;
  Z = X(check(X));  
end