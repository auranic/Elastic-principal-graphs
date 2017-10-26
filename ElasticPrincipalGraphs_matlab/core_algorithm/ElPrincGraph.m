function [NodePositions, ElasticMatrix, ReportTable]...
    = ElPrincGraph(X, NumNodes, Lambda, Mu, varargin)
% This function constructs a principal tree with NumNodes for a dataset X,
% with elasticities for edges Lamda and Elasticities for stars Mu
%
% Possible parameters
% 'InitNodePositions' - define the initial tree configuration
% 'InitElasticMatrix' - define the initial tree configuration

    np=[];
    em=[];
    verbose = 1;

    % default grammar - for a principal tree with pruning
    growGrammars = [{'bisectedge';'addnode2node'},{'bisectedge';'addnode2node'}];
    shrinkGrammars = [{'shrinkedge';'removenode'}];


    for i=1:2:length(varargin)
        if strcmpi(varargin{i},'verbose')
            verbose = varargin{i+1};
        end
        if strcmpi(varargin{i},'InitNodePositions')
            np = varargin{i+1};
        end
        if strcmpi(varargin{i},'InitElasticMatrix')
            em = varargin{i+1};
            if em~=em'
                display('ERROR: Elastic matrix must be square and symmetric');
            end
        end
        if strcmpi(varargin{i},'GrowGrammars')
            growGrammars = varargin{i+1};
        end
        if strcmpi(varargin{i},'ShrinkGrammars')
            shrinkGrammars = varargin{i+1};
        end
    end

    
if(size(em,1)==0)&(size(np,1)==0)
    edges = [1,2];
    lambdas = repmat(Lambda,size(edges,1),1);
    mus = [Mu,Mu];
    [em] = Encode2ElasticMatrix(edges,lambdas,mus);
end    
    
if(size(em,1)==0)&(size(np,1)>0)
    edges = zeros(size(np,1),2);
        for k=2:size(np,1)
            edges(k-1,1) = k-1;
            edges(k-1,2) = k;
        end
   em = Encode2ElasticMatrix(edges,Lambda*ones(size(edges,1)),Mu*ones(size(np,1)));
end


if size(np,1)==0
% default initial tree configuration - chain of nodes along the first
% principal component direction
np = zeros(size(em,1),size(X,2));
[v,u,s] = pca(X);
mn = mean(X);
v1 = v(:,1)/norm(v(:,1));
st = std(u(:,1));
delta = 2*st/(size(em,1)-1);
for k=1:size(em,1)
    np(k,:) = mn-st*v1'+delta*(k-1)*v1';
end
end

fullReport = '';

if verbose
    display(sprintf('BARCODE\tENERGY\tNNODES\tNEDGES\tNRIBS\tNSTARS\tNRAYS\tNRAYS2\tMSE MSEP\tFVE\tFVEP\tUE\tUR\tURN\tURN2\tURSD'));
end

BARCODES = {''}; 
ENERGY = zeros(NumNodes-2,1); 
NNODES = zeros(NumNodes-2,1); 
NEDGES = zeros(NumNodes-2,1); 
NRIBS = zeros(NumNodes-2,1); 
NSTARS = zeros(NumNodes-2,1); 
NRAYS = zeros(NumNodes-2,1); 
NRAYS2 = zeros(NumNodes-2,1); 
MSE = zeros(NumNodes-2,1); 
MSEP = zeros(NumNodes-2,1); 
FVE = zeros(NumNodes-2,1); 
FVEP = zeros(NumNodes-2,1); 
UE = zeros(NumNodes-2,1); 
UR = zeros(NumNodes-2,1); 
URN = zeros(NumNodes-2,1); 
URN2 = zeros(NumNodes-2,1); 
URSD = zeros(NumNodes-2,1); 

% now we grow the graph up to NumNodes
for i=1:NumNodes-2
    for k=1:size(growGrammars,2)
        [np,em] = ApplyOptimalGraphGrammarOpeation(X,np,em,growGrammars(:,k),varargin{:});
        if i==3 % this is needed to erase the star elasticity coefficient which was initially assigned to both leaf nodes,
                % one can erase this information after the number of nodes in the graph is >2
             inds = find(sum(em-diag(diag(em))>0)==1);
             
             em(inds(:),inds(:)) = 0;
        end
    end
    for k=1:size(shrinkGrammars,2)
        [np,em,ElasticEnergy] = ApplyOptimalGraphGrammarOpeation(X,np,em,shrinkGrammars(:,k),varargin{:});
    end
    
    [bcode,ENERGY(i),NNODES(i),NEDGES(i),NRIBS(i),NSTARS(i),NRAYS(i),NRAYS2(i),MSE(i),MSEP(i),FVE(i),FVEP(i),UE(i),UR(i),URN(i),URN2(i),URSD(i)] = ReportOnPrimitiveGraphEmbedment(X,np,em,varargin);
    BARCODES(i) = {bcode};
    reportString = sprintf('%s\t%f\t%i\t%i\t%i\t%i\t%i\t%i\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f',char(BARCODES(i)),ENERGY(i),NNODES(i),NEDGES(i),NRIBS(i),NSTARS(i),NRAYS(i),NRAYS2(i),MSE(i),MSEP(i),FVE(i),FVEP(i),UE(i),UR(i),URN(i),URN2(i),URSD(i));
    if verbose
        display(sprintf('%s',reportString));
    end
    
end

NodePositions = np;
ElasticMatrix = em;

    BARCODE = char(BARCODES);
    ReportTable = table(BARCODE,ENERGY,NNODES,NEDGES,NRIBS,NSTARS,NRAYS,NRAYS2,MSE,MSEP,FVE,FVEP,UE,UR,URN,URN2,URSD);
    
end

