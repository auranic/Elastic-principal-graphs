function [ElasticEnergy, MSE, EP, RP] = ...
    ComputePrimitiveGraphElasticEnergy(data, graph, part)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computes elastic energy of primitive elastic graph 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%Inputs
%   data is structure which contains all data related values. Structure
%           has following fields:
%       nPoints is number of data points.
%       dim is dimension (number of coordinates) of data
%       X is the nPoints-by-dim data matrix. Each row corresponds to one
%           data point.
%       Weights is nPoints-by-1 vector of data point weights (one weoght
%           for each point).
%       SquaredX is nPoints-by-1 vector of squared length of each data
%           vector.
%       XW = is the nPoints-by-dim matrix with data multiplied by data
%           point weights.
%   graph is structure of graph description and construction parameters.
%           Structure contains following fields:
%       nNodes is number of nodes;
%       NodePositions - is k-by-m matrix of positions of the graph nodes in
%           the same space as data. 
%       Lambda is predefined elasticity of edge
%       Mu is predefined elasticity of star
%       Lambdas is nNodes-by-nNodes symmetric matrix describing the
%           connectivity and the elastic properties of the graph.
%           Non-zero element Lambdas(i,j) is elasticity  modulo of edge
%           which connects nodes i and j. Main diagonal elements are zero.
%       Mus is vector of stars elasticities. Non-zero Mus(i) corresponds to
%       star with centre in node i.
%       TrimmingRadius is squared trimming radius, a parameter required for
%           robust principal graphs (see https://github.com/auranic/Elastic-principal-graphs/wiki/Robust-principal-graphs)
%       LocalSearch is boolean with value true corresponds to local search.
%       RadiusOfLocalSearch is graph distance to define neighbourhood for
%           local search.
%       MaxNumberOfIterations is maximal number of iteration for graph
%           embedding (usually two iterations is enough).
%       eps is minimal relative change in the node positions to be
%           considered the graph embedded (convergence criteria)
%       MaxMemorySize is the maximum memory size for the block of the
%           distance matrix when partition the data. This means that
%           maximal number of simultaneously processed data points is not
%           greater than MaxBlockSize/k, where k is number of nodes. 
%       MaxBlockSize is maximal number of simultaneously tested data
%           points: MaxBlockSize = MaxMemorySize / nNodes;
%   part is structure which contains data point associations with graph's
%           nodes. Structure contains following fields:
%       partition is nPoints-by-1 vector. partition(i) is number of node
%           which is associated with data point X(i,:). For robust version
%           of graph partition(i)==0 if point X(i,:) futher from all nodes
%           then TrimmingRadius.
%       dists is nPoints-by-1 vector. dists(i) is squared distance between
%           node with number partition(i) and data point X(i,:). For robust
%           version of graph dist(i) = TrimmingRadius^2 for all points with
%           partition(i)==0. 
%
%Outputs
%   ElasticEnergy is total elastic energy 
%   MSE is mean square error of data approximation.
%   EP is edge potential 
%   RP is harmonicity potential 

    % Calculate MSE by usage dists
    MSE = part.dists' * data.Weights / sum(data.Weights);

    % Decompose ElasticMatrix
    Lambda = triu(graph.Lambdas,1);
    StarCenterIndices = find(graph.Mus>0);
    
    % Calculate edge potential
    [row,col] = find(Lambda);
    dev = graph.NodePositions(row, :)- graph.NodePositions(col, :);
    l = Lambda(Lambda>0);
    EP = sum(l(:) .* sum(dev.^2, 2));

    % Form indeces
    indL = (graph.Lambdas) > 0;

    % Calculate harmonicity potential
    RP = 0;
    for i=1:size(StarCenterIndices,1)
        leaves = indL(:,StarCenterIndices(i));
        K = sum(leaves);
        dev = graph.NodePositions(StarCenterIndices(i),:)...
            - sum(graph.NodePositions(leaves, :)) / K;
        RP = RP + graph.Mus(StarCenterIndices(i)) * sum(dev .^ 2);
    end
    % Total energy is the sum.
    ElasticEnergy = MSE + EP + RP;
end
