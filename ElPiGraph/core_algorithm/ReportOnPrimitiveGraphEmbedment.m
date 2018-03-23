function [BARCODE, ENERGY, NNODES, NEDGES, NRIBS, NSTARS, NRAYS, NRAYS2,...
    MSE, MSEP, FVE, FVEP, UE, UR, URN, URN2, URSD] =...
    ReportOnPrimitiveGraphEmbedment(data, graph, part, computeMSEP)

% X, NodePositions, ElasticMatrix,...
%     partition, dists, computeMSEP)

% This function computes various measurements concerning a primitive graph
% embedment: 
%   BARCODE is barcode in form ...S4|S3||N, where N is number of nodes, S3
%       is number of 3-stars, S4 (S5,...) is number of four (five,...)
%       stars, etc.  
%   ENERGY is total elastic energy of graph embedment 
%           (ENERGY = MSE + UE + UR)
%   NNODES is number of nodes.
%   NEDGES is number of edges 
%   NRIBS is number of two stars (nodes with two otherr connected nodes). 
%   NSTARS is number of stars with 3 and more leaves (nodes connected
%       with central node).  
%   NRAYS2 is sum of rays minus doubled number of nodes.
%   MSE is mean square error or assessment of data approximation quality. 
%   MSEP is mean square error after piece-wise linear projection on the edges
%   FVE is fraction of explained variance. This value always between 0 
%       and 1. Greater value means higher quality of data approximation.
%   FVEP is same as FVE but computed after piece-wise linear projection
%       on the edges 
%   UE is total sum of squared edge lengths.
%   UR is total sum of star deviations from harmonicity.
%   URN is UR * nodes 
%   URN2 is UR * nodes^2
%   URSD is standard deviation of UR???
%
%Inputs:
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
%   computeMSEP is indicator of calculation of the projections not into
%       nearest node but onto nearest edge. 

    % General calculations
    TotalVariance = sum(var(data.X));
    Lambda = triu(graph.Lambdas,1);
    [row, col] = find(Lambda);
    Edges = [row, col];
    % Get barcode
    [BARCODE, N] = getPrimitiveGraphStructureBarCode(graph.Lambdas);
    % Get partition and dists if it is necessary
    if isempty(part.partition) || isempty(part.dists)
        part = PartitionDataInt(data, graph, part);
    end
    % Calculate energy
    [ENERGY, MSE, UE, UR] =...
        ComputePrimitiveGraphElasticEnergy(data, graph, part);
    NNODES = graph.nNodes;
    NEDGES = size(Edges, 1);
    nL = length(N);
    if nL > 1 
        NRIBS = N(2);
    else
        NRIBS = 0;
    end
    if nL > 2
        NSTARS = sum(N(3:end));
    else
        NSTARS = 0;
    end
    NRAYS = 0;
    NRAYS2 = 0;
    if computeMSEP
        MSEP =...
            project_point_onto_graph(data.X,  graph.NodePositions,...
            Edges, part.partition);
        FVEP = (TotalVariance - MSEP) / TotalVariance;
    else
        MSEP = NaN;
        FVEP = NaN;
    end
    FVE = (TotalVariance - MSE) / TotalVariance;
    URN = UR * NNODES;
    URN2 = UR * NNODES * NNODES;
    URSD = 0;
end

