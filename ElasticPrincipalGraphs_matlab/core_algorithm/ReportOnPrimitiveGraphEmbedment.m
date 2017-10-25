function [BARCODE,ENERGY,NNODES,NEDGES,NRIBS,NSTARS,NRAYS,NRAYS2,MSE,MSEP,FVE,FVEP,UE,UR,URN,URN2,URSD] = ReportOnPrimitiveGraphEmbedment(X,NodePositions,ElasticMatrix,varargin)
%   This function computes various measurements concerning a primitive
%   graph embedment
%   
%           BARCODE is barcode in form ...S4|S3||N, where N is number of
%               nodes, S3 is number of 3-stars, S4 (S5,...) is number of
%               four (five,...) stars, etc.
%           ENERGY is total elastic energy of graph embedment (ENERGY = MSE + UE +
%               UR)
%           NNODES is number of nodes.
%           NEDGES is number of edges 
%           NRIBS is number of two stars (nodes with two otherr connected
%               nodes). 
%           NSTARS is number of stars with 3 and more leaves (nodes
%               connected with central node). 
%           NRAYS2 is sum of rays minus doubled number of nodes.
%           MSE is mean square error or assessment of data approximation
%               quality.
%           MSEP is mean square error after piece-wise linear projection on the edges
%           FVE is fraction of explained variance. This value always
%               between 0 and 1. Greater value means higher quality of
%               data approximation.
%           FVEP is same as FVE but computed after piece-wise linear projection on the edges
%           UE is total sum of squared edge lengths.
%           UR is total sum of star deviations from harmonicity.
%           URN is UR * nodes 
%           URN2 is UR * nodes^2
%           URSD is standard deviation of UR???

    computeMSEP = 0;
    partition = [];
    
    for i=1:2:length(varargin)
        if strcmpi(varargin{i},'ComputeMSEP')
            computeMSEP = varargin{i+1};
        end
    end


Mu = diag(ElasticMatrix);
L = ElasticMatrix - diag(Mu);
Connectivities = sum(L>0);
maxc = max(Connectivities);
ee = 1:maxc;
N = histc(Connectivities,ee);
[Edges] = DecodeElasticMatrix(ElasticMatrix);
TotalVariance = sum(var(X));
BranchingFee = 0;

BARCODE = getPrimitiveGraphStructureBarCode(ElasticMatrix);

if(size(partition,1)==0)
    SquaredX = sum(X.^2, 2);
    MaxBlockSize = 100000;
    [partition,dists] = PartitionData(X,NodePositions,MaxBlockSize,SquaredX);
end


[ENERGY,MSE,UE,UR] = ComputePrimitiveGraphElasticEnergy(X,NodePositions,ElasticMatrix,partition,dists,BranchingFee);
NNODES = size(NodePositions,1);
NEDGES = size(Edges,1);
NRIBS = 0;
if size(N,2)>1 NRIBS = N(2); end;
NSTARS = 0; 
if size(N,2)>2 NSTARS = N(3); end;
NRAYS = 0;
NRAYS2 = 0;
MSEP = NaN;
if computeMSEP
    MSEP = project_point_onto_graph(X,NodePositions,Edges,partition);
end
FVE = (TotalVariance-MSE)/TotalVariance;
FVEP = NaN;
if computeMSEP
    FVEP = (TotalVariance-MSEP)/TotalVariance;
end
URN = UR*NNODES;
URN2 = UR*NNODES*NNODES;
URSD = 0;

end

