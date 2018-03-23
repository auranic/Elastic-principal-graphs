% Test of project_point_onto_graph.m
% Load data for test
load('testData.mat');
% Constants for test
nRep = 10;

%
part = PartitionData(X, np, 100000);

% test of old version
tic;
for k = 1:nRep
    [MSE1, X_projected1, EdgeIndices1, ProjectionValues1] =...
        project_point_onto_graph(X, np, ed, part);
end
t1 = toc;

% Test of Luka version
tic;
for k = 1:nRep
    [MSE2, X_projected2, EdgeIndices2, ProjectionValues2] =...
        project_point_onto_graphL(X, np, ed, part);
end
t2 = toc;

% Test of Luka version
tic;
for k = 1:nRep
    [MSE3, X_projected3, EdgeIndices3, ProjectionValues3] =...
        project_point_onto_graphL(X, np, ed, part);
end
t3 = toc;

[t1, t2, t3]