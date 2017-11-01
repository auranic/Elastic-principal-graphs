function [BARCODE, ENERGY, NNODES, NEDGES, NRIBS, NSTARS, NRAYS, NRAYS2,...
    MSE, MSEP, FVE, FVEP, UE, UR, URN, URN2, URSD] =...
    ReportOnPrimitiveGraphEmbedment(X, NodePositions, ElasticMatrix,...
    partition, dists, computeMSEP)
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
%   X is n-by-m matrix of datapints with one data point per row. n is
%       number of data points and m is dimension of data space.
%   NodePositions is k-by-m matrix of embedded coordinates of graph nodes,
%       where k is number of nodes and m is dimension of data space.
%   ElasticMatrix is k-by-k matrix of nodes connectivity: 
%       ElsticMatrix(i,i) > 0 if node i is centre of star and zero otherwise
%       ElsticMatrix(i,j) > 0 if there is edge between nodes i and j. In
%       this case ElsticMatrix(i,j) is elasticity modulo of edge from i to j.
%   partition is n-by-1 vector. partition(i) is number of node which is
%       associated with data point X(i,:).
%   ComputeMSEP is non-zero to compute MSEP and zero otherwise.
%

    % General calculations
    TotalVariance = sum(var(X));
    BranchingFee = 0;
    Lambda = triu(ElasticMatrix,1);
    [row, col] = find(Lambda);
    Edges = [row, col];
    % Get barcode
    [BARCODE, N] = getPrimitiveGraphStructureBarCode(ElasticMatrix);
    % Get partition and dists if it is necessary
    if(size(partition,1)==0 || size(dists,1)==0)
        [partition, dists] = PartitionData(X, NodePositions, 100000,...
            sum(X.^2, 2));
    end
    % Calculate energy
    [ENERGY, MSE, UE, UR] =...
        ComputePrimitiveGraphElasticEnergy(X, NodePositions,...
        ElasticMatrix, partition, dists, BranchingFee);
    NNODES = size(NodePositions, 1);
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
        MSEP = project_point_onto_graph(X,  NodePositions,  Edges, partition);
        FVEP = (TotalVariance-MSEP)/TotalVariance;
    else
        MSEP = NaN;
        FVEP = NaN;
    end
    FVE = (TotalVariance-MSE)/TotalVariance;
    URN = UR * NNODES;
    URN2 = UR * NNODES * NNODES;
    URSD = 0;
end

