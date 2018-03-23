function [NodePositions, ElasticMatrix, ReportTable]...
    = ElPrincGraph(X, NumNodes, Lambda, Mu, varargin)
% This function constructs a principal tree with NumNodes for a dataset X,
% with elasticities for edges Lamda and Elasticities for stars Mu
%
% Outputs
%   np is positions of empbedded nodes.
%   em is elastic matrix of the best selected new graph.
%   ReportTable is structure with 16 fields:
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
%       MSEP is mean square error for projections to edge.
%       FVE is fraction of explained variance. This value always between 0
%           and 1. Greater value means higher quality of data approximation.
%       FVEP is fraction variance unexplained for projections to edge.
%       UE is elastic energy for edges stretching.
%       UR is elastic energy of deviation from harmonicity.
%       URN is UR * nodes 
%       URN2 is UR * nodes^2
%       URSD is standard deviation of UR
%
%Inputs
%   X - is the n-by-m data matrix. Each row corresponds to one data point.
%   NumNodes maximal number of nodes in graph.
%   Lambda is edge elasticity potential
%   Mu is star elasticity potential
%   varargin contains Name, Value pairs. Names can be: 
%   'InitNodePositions' is k-by-m matrix of the initial tree node positions.
%   'InitElasticMatrix' is initial elastic matrix: k-by-k symmetric matrix
%       describing the connectivity and the elastic properties of the
%       graph. Star elasticities (mu coefficients) are presented on the
%       main diagonal (non-zero entries only for star centers), and the
%       edge elasticity moduli are presented out of diagonal. 
%   'verbose' with 1/0 is to display/hide the energy values at each
%       iteration and in the end of the process. 
%   'PointWeights' with n-by-1 vector of data point weights
%   'GrowGrammars' is cell array of string which contain grammar names:
%       'addnode2node'      adds a node to each graph node
%       'bisectedge'        adds nodt to the middle of each edge
%   'ShrinkGrammars' is cell array of string which contain grammar names:
%       'removenode'        removes terminal node
%       'shrinkedge'        removes edge and glues two ends of this edge
%   'ComputeMSEP' with 1/0 is indicator of calculation of the projections
%       not into nearest node but onto nearest edge.
%   'TrimmingRadius' is trimming radius, a parameter required for robust
%       principal graphs (see https://github.com/auranic/Elastic-principal-graphs/wiki/Robust-principal-graphs)
%   'MaxMmorySize' with integer number which is maximum memory size for the
%       block of the distance matrix when partition the data. This means
%       that maximal number of simultaneously processed data points is not
%       greater than MaxBlockSize/k, where k is number of nodes.
%   'LocalSearch' specifies local test for optimizing only a subset of the
%       nodes. Integer number defines radius of neghbourhood to use.
%   'MaxNumberOfIterations' with integer number which is maximum number of
%       iterations for EM algorithm. 
%   'eps' with real number which is minimal relative change in the node
%       positions to be considered the graph embedded (convergence criteria) 

    % Prepare structures for data
    data.X = X;
    [data.nPoints, data.dim] = size(X);
    data.Weights = [];
    data.SquaredX = sum(X .^ 2, 2);
    data.XW = []; % Will be completed after the weights identification.

    % Prepare structures for graph
    graph.nNodes = 0;
    graph.NodePositions = [];
    graph.Lambda = Lambda;
    graph.Mu = Mu;
    graph.Lambdas = [];
    graph.Mus = [];
    graph.TrimmingRadius = Inf;     % We remember squared trimming radius!
    graph.LocalSearch = false;
    graph.RadiusOfLocalSearch = 0;
    graph.NodesSubSet = [];
    graph.MaxNumberOfIterations = 10;
    graph.eps = 0.01;
    graph.MaxBlockSize = 0;          % MaxMemorySize/nNodes;
    graph.MaxMemorySize = 10000000;  % 80M is maximal size to use for 
                                     % distance matrix calculation 

    % Prepare structures for partition
    part.partition = [];
    part.dists = [];
    
    % Prepare variables for local parameters
    verbose = true;
    computeMSEP = 0;

    % default grammar - for a principal tree with pruning
    growGrammars = [{'bisectedge';'addnode2node'},{'bisectedge';'addnode2node'}];
    shrinkGrammars = {'shrinkedge';'removenode'};
    
    graph.PenalizedEnergy = false;

    % Parse input arguments
    for i=1:2:length(varargin)
        if strcmpi(varargin{i},'verbose')
            verbose = varargin{i + 1};
        elseif strcmpi(varargin{i},'ComputeMSEP')
            ComputeMSEP = varargin{i + 1};
        elseif strcmpi(varargin{i},'GrowGrammars')
            growGrammars = varargin{i + 1};
        elseif strcmpi(varargin{i},'ShrinkGrammars')
            shrinkGrammars = varargin{i + 1};
        elseif strcmpi(varargin{i},'InitNodePositions')
            graph.NodePositions = varargin{i + 1};
        elseif strcmpi(varargin{i},'InitElasticMatrix')
            graph.Lambdas = varargin{i + 1};
            if graph.Lambdas ~= graph.Lambdas' %#ok<BDSCA>
                error('ERROR: Elastic matrix must be square and symmetric');
            end
            graph.Mus = diag(graph.Lambdas);
            graph.Lambdas = graph.Lambdas - diag(graph.Mus);
        elseif strcmpi(varargin{i},'LocalSearch')
            graph.LocalSearch = 1;
            graph.RadiusOfLocalSearch = varargin{i + 1};
        elseif strcmpi(varargin{i},'TrimmingRadius')
            graph.TrimmingRadius = varargin{i + 1};
            if isdeployed
                graph.TrimmingRadius = str2double(graph.TrimmingRadius);
            end
            graph.TrimmingRadius = graph.TrimmingRadius .^ 2;
        elseif strcmpi(varargin{i},'MaxNumberOfIterations')
            graph.MaxNumberOfIterations = varargin{i + 1};
            if isdeployed
                graph.MaxNumberOfIterations =...
                    str2double(graph.MaxNumberOfIterations);
            end
        elseif strcmpi(varargin{i},'eps')
            graph.eps = varargin{i + 1};
            if isdeployed
                graph.eps = str2double(graph.eps);
            end
        elseif strcmpi(varargin{i},'MaxMemorySize')
            graph.MaxMemorySize = varargin{i + 1};
            if isdeployed
                graph.MaxMemorySize = str2double(MaxBlockSize);
            end
        elseif strcmpi(varargin{i},'PointWeights')
            data.PointWeights = varargin{i + 1};
        elseif strcmpi(varargin{i},'BranchingControls')
            graph.PenalizedEnergy = true;
            graph.BranchingControls = varargin{i + 1};
        end
        
    end

    % Check the necessity of report
    if nargout<3 && ~verbose
        report = false;
    else
        report = true;
    end
    
    % Check PointWeights
    if isempty(data.Weights)
        data.Weights = ones(data.nPoints, 1);
        data.XW = data.X;
    else
        if length(data.PointWeights) ~= data.nPoints
            error('Number of point weights must be the same as nmber of points');
        end
        data.PointWeights = data.PointWeights(:);
        data.XW = bsxfun(@times,X, data.PointWeights);
    end
    
    % Check consistency of initial elastic matrix and node positions
    if ~isempty(graph.NodePositions) && ~isempty(graph.Lambdas)
        if size(graph.NodePositions,1) ~= size(graph.Lambdas,1)
            error(['Number of nodes (number of rows in NodePositions)'...
                ' must coinside with dimension of Elastic matrix']);
        end
    end
    
    % Fix number of nodes as it is known at the moment
    graph.nNodes = size(graph.NodePositions, 1);
    % Create initial ElasticMatrices if it is not specified
    if isempty(graph.Lambdas)
        % If initial node position is not specified then create place for
        % two nodes
        if graph.nNodes == 0
            edges = [1,2];
            graph.nNodes = 2;
        else
            % Otherwise connect each pair of consequent nodes
            edges = [(1:graph.nNodes-1)', (2:graph.nNodes)'];
        end
        % Create initial ElasticMatrices for specified structure of graph
        graph = MakeUniformElasticMatrix(edges, graph);
    end
    
    % Check if the NodePositions are defined 
    if isempty(graph.NodePositions)
        % default initial tree configuration - chain of nodes along the
        % first principal component direction
        % Compute mean of X
        mn = mean(X);
        % Compute the first PC
        XM = bsxfun(@minus, X, mn);
        [~,~,v] = svds(XM, 1);
        % Compute standard projections onto the first principal component
        % and standard deviation of it.
        st = std(XM * v);
        % Dispersion will be 2 standard deviation. This means that step is
        delta = 2 * st / (graph.nNodes - 1);
        graph.NodePositions = bsxfun(@plus, mn - st * v',...
            (delta * (0:graph.nNodes - 1)') * v');
        % Clear memory
        clear mn XM v st delta  
    end
    
    %Now we can recalculate MaxBlockSize
    graph.MaxBlockSize = floor(graph.MaxMemorySize / graph.nNodes);
    
    if verbose
        display(sprintf(['BARCODE\tENERGY\tNNODES\tNEDGES\tNRIBS\tNSTARS'...
            '\tNRAYS\tNRAYS2\tMSE MSEP\tFVE\tFVEP\tUE\tUR\tURN\tURN2\tURSD']));
    end
    
    % Initialise arrays for report
    if report
        BARCODES = {''};
        ENERGY = zeros(NumNodes - graph.nNodes, 1); 
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

    % First, we optimize the graph without any growth
    
    %[graph, part, ~] = PrimitiveElasticGraphEmbedment(data, graph, part);
    [graph, part] = ApplyOptimalGraphGrammarOperation(data, graph, part,{});

    
    % Now we grow the graph up to NumNodes
    i=1; % Number of row of current graph in report
    while (graph.nNodes < NumNodes)
        % Growing
        for k=1:size(growGrammars,2)
            [graph, part] =...
                ApplyOptimalGraphGrammarOperation(data, graph, part,...
                growGrammars(:, k));
        end
        % Shrinking
        for k=1:size(shrinkGrammars, 2)
            [graph, part] =...
                ApplyOptimalGraphGrammarOperation(data, graph, part,...
                shrinkGrammars(:,k));
        end
        % form report if it is necessary
        if report
            % Compute data for report
            [bcode, ENERGY(i), NNODES(i), NEDGES(i), NRIBS(i), NSTARS(i),...
                NRAYS(i), NRAYS2(i), MSE(i), MSEP(i), FVE(i), FVEP(i),...
                UE(i), UR(i), URN(i), URN2(i), URSD(i)] =...
                ReportOnPrimitiveGraphEmbedment(data, graph, part,...
                computeMSEP);
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
        i = i - 1;
        BARCODE = char(BARCODES(1:i));
        ENERGY = ENERGY(1:i);
        NNODES = NNODES(1:i);
        NEDGES = NEDGES(1:i); 
        NRIBS = NRIBS(1:i); 
        NSTARS = NSTARS(1:i); 
        NRAYS = NRAYS(1:i); 
        NRAYS2 = NRAYS2(1:i); 
        MSE = MSE(1:i); 
        MSEP = MSEP(1:i); 
        FVE = FVE(1:i); 
        FVEP = FVEP(1:i); 
        UE = UE(1:i); 
        UR = UR(1:i); 
        URN = URN(1:i); 
        URN2 = URN2(1:i); 
        URSD = URSD(1:i);
        ReportTable = struct('BARCODE', BARCODE, 'ENERGY', ENERGY,...
            'NNODES', NNODES, 'NEDGES', NEDGES, 'NRIBS', NRIBS,...
            'NSTARS',NSTARS, 'NRAYS', NRAYS, 'NRAYS2', NRAYS2,...
            'MSE', MSE, 'MSEP', MSEP, 'FVE', FVE, 'FVEP', FVEP,...
            'UE', UE, 'UR', UR, 'URN', URN, 'URN2', URN2, 'URSD', URSD);
    end
    NodePositions = graph.NodePositions;
    ElasticMatrix = graph.Lambdas + diag(graph.Mus);
end

function graph = MakeUniformElasticMatrix(Edges, graph)
% Creates an elastic matrix for a primitiveGraph defined by its vector of
% edges
% The elasticities of the edges are equal to Lambda for all edges
% and to Mu for all stars.
%
% Inputs:
%   Edges is n-by-2 array with number of one edge's node in column 1 and
%       number of another node of the same edge in the column 2.
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

    % Create ElasticMatrix
    graph.Lambdas = zeros(graph.nNodes, graph.nNodes);
    % Transform matrix of Edges into index of ElasticMatrix
    ind = sub2ind(size(graph.Lambdas), Edges(:,1), Edges(:,2));
    % Fill all edge elasticities
    graph.Lambdas(ind) = graph.Lambda;
    graph.Lambdas = graph.Lambdas + graph.Lambdas';
    % Calculate connectivities
    Connect = sum(graph.Lambdas > 0);
    % Identify stars
    ind  = Connect > 1;
    % Generate vector od star centres
    graph.Mus = zeros(graph.nNodes, 1);
    % Specify star elasticities
    graph.Mus(ind) = graph.Mu;
end