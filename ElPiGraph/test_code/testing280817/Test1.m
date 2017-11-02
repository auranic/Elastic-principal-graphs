%Test of work
%General constants
nData = 1000000;
dim = 100;
nNodes = 50;
X = rand(nData,dim);
% Randomly select nNodes nodes
ind = randsample(nData,nNodes);
NodePositions = X(ind,:);

XSquared = sum(X.^2,2);

tic;
[partition,dists] = PartitionData(X,NodePositions,100000,XSquared);
toc