%Test of two methods
nRep = 100;
nRows = 10000;
nCols = 50;
nNodes = 5;
% Generate data
data = rand(nRows, nCols);
% Partitions
part = randi(nNodes, nRows, 1);
% Weights
weights = rand(nRows, 1);
%Test of the null method
tic;
for k = 1:nRep
    [cent, rad] = ComputeWeightedAverage(data, part, weights, nNodes);
end
toc
%Test of the first method
tic;
for k = 1:nRep
    [cent, rad] = ComputeWeightedAverage4(data, part, weights, nNodes);
end
toc
%Test of the second method
tic;
for k = 1:nRep
    [cent2, rad2] = ComputeWeightedAverage2(data, part, weights, nNodes);
end
toc
%Test of the third method
tic;
for k = 1:nRep
    [cent2, rad2] = ComputeWeightedAverage3(data, part, weights, nNodes);
end
toc