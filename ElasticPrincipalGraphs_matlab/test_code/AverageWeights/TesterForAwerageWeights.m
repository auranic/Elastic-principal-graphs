%Tester for Awerage weights
%Old version
nRep = 100;
tic;
for k = 1:nRep
    [NodeClusterCenters, NodeClusterRelativeSize] =...
        ComputeWeightedAverage(X, partition, PointWeights, NumberOfNodes);    
end
r1 = toc/nRep;
%New version
tic;
for k = 1:nRep
    [NodeClusterCenters1, NodeClusterRelativeSize1] =...
        ComputeWeightedAverageM(X, partition, PointWeights, NumberOfNodes);    
end
r2 = toc/nRep;
%New version
tic;
for k = 1:nRep
    [NodeClusterCenters2, NodeClusterRelativeSize2] =...
        ComputeWeightedAverageM(X, partition, PointWeights, NumberOfNodes);    
end
r3 = toc/nRep;
[r1 r2 r3 isequal(NodeClusterCenters,NodeClusterCenters1,NodeClusterCenters2) isequal(NodeClusterRelativeSize,NodeClusterRelativeSize1,NodeClusterRelativeSize2)]
