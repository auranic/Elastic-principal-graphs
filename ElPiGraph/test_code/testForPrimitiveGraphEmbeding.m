X = load('./test_data/tree23/tree23.data');
%X = load('./test_data/mosaic/mosaic.txt');

% inflate the number of points
inflationFactor = 100;
X1 = zeros(size(X,1)*inflationFactor,size(X,2));
k=1;
STDV = std(X);
if 0
for i=1:size(X,1)
    for j=1:inflationFactor
        r = rand(1,size(X,2));
        p = 0.5*STDV.*r;
        X1(k,:) = X(i,:)+p;
        k=k+1;
    end
end
X = X1;
end

Npoints = size(X,1)

%Edges = [1,4;2,4;3,4];

Edges = [1,2;2,3;3,4;4,5;5,6;6,7;7,8;4,9;9,10;10,11;11,12];
%Edges = [1,2;2,3;3,4;4,5;5,6;6,7;7,8];

NNodes = max(max(Edges));

ElasticMatrix = MakeUniformElasticMatrix(Edges,0.01,0.1);
%ElasticMatrix = MakeUniformElasticMatrix(Edges,0,0);

%NodePositions = X(randperm(Npoints,NNodes),:);

NodePositions = zeros(size(ElasticMatrix,1),size(X,2));
[v,u,s] = pca(X);
mn = mean(X);
v1 = v(:,1)/norm(v(:,1))
%v1 = [1;0;0]
st = std(u(:,1));
delta = 2*st/(size(ElasticMatrix,1)-1);
for k=1:size(ElasticMatrix,1)
    NodePositions(k,:) = mn-st*v1'+delta*(k-1)*v1';
end



tic
%[EmbeddedNodePositions,~,partition] = PrimitiveElasticGraphEmbedment(X, NodePositions, ElasticMatrix,'MaxNumberOfIterations',100,'verbose',0);

% testing local optimization
Xsquared = sum(X.^2,2); 
[partition] = PartitionData(X,NodePositions,100000,Xsquared);
PointWeights = ones(size(X,1),1);
NodeSubSet = [1:8];
%NodeSubSet = [3];
[NodeClusterCenters, NodeClusterRelativeSize] = ComputeWeightedAverage(X, partition, PointWeights, size(NodePositions,1));
[EmbeddedNodePositions,~,partition] = PrimitiveElasticGraphEmbedment(X, NodePositions, ElasticMatrix,'MaxNumberOfIterations',100,...
   'verbose',1,...
   'Local',struct('Nodes',NodeSubSet,'Partition',partition));

%[EmbeddedNodePositions,~,partition] = PrimitiveElasticGraphEmbedment(X, NodePositions, ElasticMatrix,'MaxNumberOfIterations',100,'verbose',1);

toc

plot(X(:,1),X(:,2),'k.'); hold on; drawGraph2D(EmbeddedNodePositions,Edges);

