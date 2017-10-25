%Test of work
%General constants
nData = 100000;
dim = 100;
nNodes = 50;
Radius = Inf;
threshold = 1.4;
% Generate random data
X = rand(nData,dim);
% Randomly select nNodes nodes
ind = randsample(nData,nNodes);
NodePositions = X(ind,:);
% Create elastic matrix
ElasticMatrix = rand(nNodes);
%Transfer to symmetry
ElasticMatrix = ElasticMatrix + ElasticMatrix';
%Kill small values
ElasticMatrix(ElasticMatrix<threshold) = 0;
%Kill diagonel
ElasticMatrix = ElasticMatrix - diag(diag(ElasticMatrix));
%Create diagonal
ElasticMatrix = ElasticMatrix + diag(sum(ElasticMatrix>0)>1);
%Test function
tic;
[EmbeddedNodePositions, ElasticEnergy, partition, MSE,EP,RP] =...
    PrimitiveElasticGraphEmbedment(X, NodePositions, ElasticMatrix);
toc