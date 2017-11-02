function [np, em, ReportTable]...
    = ElPrincGraph(X, NumNodes, Lambda, Mu, varargin)
% This function constructs a principal tree with NumNodes for a dataset X,
% with elasticities for edges Lamda and Elasticities for stars Mu
%
% Outputs
%   np is positions of empbedded nodes.
%   em is elastic matrix of the best selected new graph.
%   ReportTable is table with 17 columns:
%       STEP is number of iteration
%       BARCODE is barcode in form ...S4|S3||N, where N is number of
%           nodes, S3 is number of three stars, S4 (5,...) is number of
%           four (five,...) stars.
%       ENERGY is total energy of graph embedment (ENERGY = MSE + UE + UR)
%       NNODES is number of nodes.
%       NEDGES is number of edges (NEDGES = NNODES - 1).
%       NRIBS is number of two stars (nodes with two otherr connected nodes). 
%       NSTARS is number of stars with 3 and more leaves (nodes connected
%           with central node).  
%       NRAYS2 is sum of rays minus doubled number of nodes.
%       MSE is mean square error or assessment of data approximation
%           quality. 
%       MSEP is project to edge ???
%       FVE is fraction of explained variance. This value always between 0
%           and 1. Greater value means higher quality of data approximation.
%       FVEP is ???
%       UE is elastic energy for edges stretching.
%       UR is elastic energy of deviation from harmonicity.
%       URN is UR * nodes 
%       URN2 is UR * nodes^2
%       URSD is standard deviation of UR???
%
%Inputs
%   X - is the n-by-m data matrix. Each row corresponds to one data point.
%   NumNodes maximal number of nodes in graph.
%   Lambda is edge elasticity potential
%   Mu is star elasticity potential
%   varargin contains Name, Value pairs. Names can be: 
%       'InitNodePositions' is k-by-m matrix of the initial tree node
%           positions.
%       'InitElasticMatrix' is initial elastic matrix: k-by-k symmetric
%           matrix describing the connectivity and the elastic properties
%           of the graph. Star elasticities (mu coefficients) are presented
%           on the main diagonal (non-zero entries only for star centers),
%           and the edge elasticity moduli are presented out of diagonal.
%       'verbose' with 1/0 is to display/hide the energy values at each
%           iteration and in the end of the process.
%       Other Name/value pairs can be presented to transfer to called
%       functions like ApplyOptimalGraphGrammarOpeation.

    % Parse 
    np=[];
    em=[];
    verbose = 1;
    ComputeMSEP = 0;

    % default grammar - for a principal tree with pruning
    growGrammars = [{'bisectedge';'addnode2node'},{'bisectedge';'addnode2node'}];
    shrinkGrammars = {'shrinkedge';'removenode'};


    for i=1:2:length(varargin)
        if strcmpi(varargin{i},'verbose')
            verbose = varargin{i+1};
        elseif strcmpi(varargin{i},'InitNodePositions')
            np = varargin{i+1};
        elseif strcmpi(varargin{i},'InitElasticMatrix')
            em = varargin{i+1};
            if em ~= em' %#ok<BDSCA>
                error('ERROR: Elastic matrix must be square and symmetric');
            end
        elseif strcmpi(varargin{i},'GrowGrammars')
            growGrammars = varargin{i+1};
        elseif strcmpi(varargin{i},'ShrinkGrammars')
            shrinkGrammars = varargin{i+1};
        elseif strcmpi(varargin{i},'ComputeMSEP')
            ComputeMSEP = varargin{i+1};
        end
    end

    % Check the necessity of report
    if nargout<3 && ~verbose
        report = false;
    else
        report = true;
    end
    
    
    CurrentNumberOfNodes = size(np, 1);
    % Create initial ElasticMatrices if it is not specified
    if size(em,1) == 0
        % If initial node position is not specified then create place for
        % two nodes
        if CurrentNumberOfNodes == 0
            edges = [1,2];
        else
            edges = [(1:CurrentNumberOfNodes-1)', (2:CurrentNumberOfNodes)];
        end
        % Create initial ElasticMatrices for specified structure of graph
        em = MakeUniformElasticMatrix(edges, Lambda, Mu);
    end

    if size(np, 1) == 0
        % default initial tree configuration - chain of nodes along the
        % first principal component direction
        CurrentNumberOfNodes = size(em, 1);
        % Compute the first PC
        [~,~,v] = svds(X,1);
        % Compute mean of X
        mn = mean(X);
        % Compute standard projections onto the first principal component
        % and standard deviation of it.
        st = std(X * v);
        % Dispersion will be 2 standard deviation. This means that step is
        delta = 2 * st / (CurrentNumberOfNodes - 1);
        np = bsxfun(@plus, mn - st * v',...
            (delta * (0:CurrentNumberOfNodes - 1)') * v');
    end
    
    % Just in case redefine this value
    CurrentNumberOfNodes = size(np, 1);

    %Control the existance of Mu in initial Elastic matrix em
    UR = diag(em);
    if sum(UR>0) == 0
        % There are no non-zero elements on main diagonal. Create it
        em = em + diag(Mu*ones(CurrentNumberOfNodes, 1));
    end

    if verbose
        display(sprintf(['BARCODE\tENERGY\tNNODES\tNEDGES\tNRIBS\tNSTARS'...
            '\tNRAYS\tNRAYS2\tMSE MSEP\tFVE\tFVEP\tUE\tUR\tURN\tURN2\tURSD']));
    end
    % Initialise arrays for report
    if report
        BARCODES = {''};
        ENERGY = [];%zeros(NumNodes - CurrentNumberOfNodes, 1); 
        NNODES = ENERGY; 
        NEDGES = ENERGY; 
        NRIBS = ENERGY; 
        NSTARS = ENERGY; 
        NRAYS = ENERGY; 
        NRAYS2 = ENERGY; 
        MSE = ENERGY; 
        MSEP = ENERGY; 
        FVE = ENERGY; 
        FVEP = ENERGY; 
        UE = ENERGY; 
        UR = ENERGY; 
        URN = ENERGY; 
        URN2 = ENERGY; 
        URSD = ENERGY;
    end

    % Now we grow the graph up to NumNodes
    i=1; % Number of row of current graph in report
    while(size(np,1)<NumNodes)
        % Growing
        for k=1:size(growGrammars,2)
            [np, em, partition, dists] =...
                ApplyOptimalGraphGrammarOpeation(X, np, em,...
                growGrammars(:, k), varargin{:});
        end
        % Shrinking
        for k=1:size(shrinkGrammars, 2)
            [np, em, partition, dists] = ...
                ApplyOptimalGraphGrammarOpeation(X, np, em,...
                shrinkGrammars(:,k), varargin{:});
        end
        % form report if it is necessary
        if report
            % Compute data for report
            [bcode, ENERGY(i,1), NNODES(i,1), NEDGES(i,1), NRIBS(i,1), NSTARS(i,1),...
                NRAYS(i,1), NRAYS2(i,1), MSE(i,1), MSEP(i,1), FVE(i,1), FVEP(i,1),...
                UE(i,1), UR(i,1), URN(i,1), URN2(i,1), URSD(i,1)] =...
                ReportOnPrimitiveGraphEmbedment(X, np, em, partition,...
                dists, ComputeMSEP);
            BARCODES(i) = {bcode};
            if verbose
                display(sprintf('%s\t%f\t%i\t%i\t%i\t%i\t%i\t%i\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f',...
                    char(BARCODES(i)), ENERGY(i), NNODES(i), NEDGES(i),...
                    NRIBS(i), NSTARS(i), NRAYS(i), NRAYS2(i), MSE(i),...
                    MSEP(i), FVE(i), FVEP(i), UE(i), UR(i), URN(i),...
                    URN2(i), URSD(i)));
            end
        end    
        i=i+1;
    end
    % Form report if it is necessary
    if nargout > 2
        BARCODE = char(BARCODES);
        ReportTable = table(BARCODE, ENERGY, NNODES, NEDGES, NRIBS, NSTARS,...
            NRAYS, NRAYS2, MSE, MSEP, FVE, FVEP, UE, UR, URN, URN2, URSD);
    end
end
