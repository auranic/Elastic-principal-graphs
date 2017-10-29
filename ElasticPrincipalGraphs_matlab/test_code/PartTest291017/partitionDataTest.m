%Partition data test
%Required arrays.
%X is array of datapoints.
%np is array of nodes' positions.
    MaxBlockSize = 100000;
    SqX = sum(X.^2,2);
    TrimmingRadius = Inf;
    
    nReps = 100;
    tic;
    for K = 1:nReps
        [partition1, dists1] = ...
            PartitionData(X, np, MaxBlockSize, SqX, TrimmingRadius);
    end
    t1 = toc / nReps;
        tic;
    for K = 1:nReps
        [partition2] = ...
            PartitionData(X, np, MaxBlockSize, SqX, TrimmingRadius);
    end
    t2 = toc / nReps;

    nReps = 100;
    tic;
    for K = 1:nReps
        [partition3, dists3] = ...
            PartitionDataM(X, np, MaxBlockSize, SqX, TrimmingRadius);
    end
    t3 = toc / nReps;
        tic;
    for K = 1:nReps
        [partition4] = ...
            PartitionDataM(X, np, MaxBlockSize, SqX, TrimmingRadius);
    end
    t4 = toc / nReps;
    
    [t1, t2, t3, t4, isequal(partition1,partition2,partition3,partition4), isequal(dists1,dists3)]