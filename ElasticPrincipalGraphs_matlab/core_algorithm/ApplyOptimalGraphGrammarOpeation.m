function [NodePositions2,ElasticMatrix2, ElasticEnergy] = ApplyOptimalGraphGrammarOpeation(X,NodePositions,ElasticMatrix,operationtypes,varargin)
% this functinon applies the most optimal graph grammar operation of operationtype
% the embedment of an elastic graph described by ElasticMatrix

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


k=1;
%tic
for i=1:size(operationtypes,1)
    %display(sprintf('Operation type = %s',char(operationtypes(i))));
    [NodePositionArray, ElasticMatrices, NodeIndicesArray] = GraphGrammarOperation(NodePositions,ElasticMatrix,X,char(operationtypes(i)));
    for j=1:size(NodePositionArray,3)
        NodePositionArrayAll(:,:,k) = NodePositionArray(:,:,j);
        ElasticMatricesAll(:,:,k) = ElasticMatrices(:,:,j);
        NodeIndicesArrayAll(:,k) = NodeIndicesArray(:,j);
        k = k+1;
    end
    %display(sprintf('Time spent for applying %s: %f',char(operationtypes(i)),toc));
end


minEnergy = realmax;
k = -1;

% we compute these things here in order not to recompute for each graph embedment
SquaredX = sum(X.^2, 2);


LocalInfo = struct();
LocalInfo.Partition = [];

if LocalSearch
    [LocalInfo.Partition] = PartitionData(X, NodePositions, MaxBlockSize, SquaredX, TrimmingRadius);
end


%tic
for i=1:size(ElasticMatricesAll,3)
    if LocalSearch
        NodeIndices = NodeIndicesArrayAll(:,i);
        sd = fast_setdiff1(1:size(NodePositionArrayAll(:,:,i),1),NodeIndices);
        if(size(sd,2)==0)
            sd = [size(NodePositions,1)];
            sd1 = GetNeighbourhoodOnTheGraph(ElasticMatrix,sd,1);
            sd = fast_setdiff1(sd1,sd);
        end
        LocalInfo.NodeSubSet = GetNeighbourhoodOnTheGraph(ElasticMatricesAll(:,:,i),sd,RadiusOfLocalSearch);
        [np,ElasticEnergy] = PrimitiveElasticGraphEmbedment(X, NodePositionArrayAll(:,:,i), ElasticMatricesAll(:,:,i),'verbose',0,'SquaredX',SquaredX,'Local',LocalInfo);
    else
        [np,ElasticEnergy] = PrimitiveElasticGraphEmbedment(X, NodePositionArrayAll(:,:,i), ElasticMatricesAll(:,:,i),'verbose',0,'SquaredX',SquaredX);
    end
    nps(:,:,i) = np(:,:);
    if ElasticEnergy<minEnergy
        minEnergy = ElasticEnergy;
        k = i; 
    end
end

if LocalSearch
    [NodePositions2,ElasticEnergy] = PrimitiveElasticGraphEmbedment(X, NodePositionArrayAll(:,:,k), ElasticMatricesAll(:,:,k),'verbose',0,'SquaredX',SquaredX);
else
    NodePositions2(:,:) = nps(:,:,k);
    ElasticEnergy = minEnergy;
end
    ElasticMatrix2 = ElasticMatricesAll(:,:,k);

%display(sprintf('Time spent for optimizing graphs: %f',toc));


end

function Z = fast_setdiff1(X,Y)
if ~isempty(X)&&~isempty(Y)
  X = X+1;
  Y = Y+1;
  check = false(1, max(max(X), max(Y)));
  check(X) = true;
  check(Y) = false;
  Z = X(check(X));  
  Z = Z-1;
else
  Z = X;
end
end

