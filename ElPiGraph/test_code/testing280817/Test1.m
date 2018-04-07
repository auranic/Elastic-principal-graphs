%Test of work
%General constants
nData = 1000000;
dim = 100;
nNodes = 50;
X = rand(nData,dim);
% Randomly select nNodes nodes
ind = randsample(nData,nNodes);
NodePositions = X(ind,:);

tic;
XSquared = sum(X.^2,2);
toc;

tic;
[partition,dists] = PartitionData(X,NodePositions,100000,XSquared);
toc

tic;
[idx,dist] = knnsearch(NodePositions,X,'k',1);
toc;