%Partition data test
%Required arrays.
%X is array of datapoints.
%np is array of nodes' positions.

    np = np(1:10,:);

    data.X = X;
    [data.nPoints, data.dim] = size(X);
    data.Weights = ones(data.nPoints,1);
    data.SquaredX = sum(X .^ 2, 2);
    data.XW = X; 

    graph.nNodes = size(np, 1);
    graph.NodePositions = np;
    graph.Lambda = 0.01;
    graph.Mu = 0.1;
    graph.Lambdas = zeros(graph.nNodes);
    graph.Mus = zeros(graph.nNodes,1);
    graph.TrimmingRadius = Inf;
    graph.LocalSearch = false;
    graph.RadiusOfLocalSearch = 0;
    graph.NodesSubSet = [];
    graph.MaxNumberOfIterations = 10;
    graph.eps = 0.01;
    graph.MaxMemorySize = 10000000;  % 80M is maximal size to use for 
                                     % distance matrix calculation 
    graph.MaxBlockSize = graph.MaxMemorySize / graph.nNodes;

    part.partition = zeros(data.nPoints, 1);
    part.dists = zeros(data.nPoints, 1);
    
    bigTest = 100;
    t = zeros(bigTest, 3);
    for b = 1:bigTest
    
    nReps = 100;
    tic;
    for K = 1:nReps
        [partition1, dists1] = ...
            PartitionData(X, np, graph.MaxBlockSize, data.SquaredX, graph.TrimmingRadius);
    end
    t(b, 1) = toc;% / nReps;
        tic;
    for K = 1:nReps
        part = ...
            PartitionDataInt(data, graph, part);
    end
    t(b, 2) = toc;% / nReps;
        tic;
    for K = 1:nReps
        part2 = ...
            PartitionDataInt2(data, graph, part);
    end
    t(b, 3) = toc;% / nReps;
    b
    end
    