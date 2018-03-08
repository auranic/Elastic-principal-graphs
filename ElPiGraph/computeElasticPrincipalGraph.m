function [NodePositions, Edges, ReportTable] =...
    computeElasticPrincipalGraph(data,NumNodes,varargin)
%computeElasticPrincipalGraph calculate elastic principal graph. 
%If grammar parameters 'GrowGrammar' and 'ShrinkGrammar' are not specified
%then by default the algorithm constructs a principal tree with default
%elasticity coefficients (Lambda=0.01,Mu=0.1) 
%
%
%Theory and some examples of computeElasticPrincipalTree usage can be found
%in https://github.com/auranic/Elastic-principal-graphs/wiki/Basic-use-of-Elastic-Principal-Graphs-Matlab-package
%
%Usage
%   [NodePositions, Edges] = computeElasticPrincipalGraph(data, NumNodes,
%   ParameterSet) returns 
%       NodePositions is k-by-m matrix of coordinates of k nodes in m
%           dimensional space for n-by-m matrix data (each row of matrix
%           contains one observation) and at most NumNodes nodes k.
%       Edges is k-by-2 matrix of integers. Edges(i,1) and Edges(i,2)
%           specify numbers of two vertex of i-th edge.
%   Input features:
%       data is n-by-m matrix data (each row of matrix contains one
%           observation).
%       NumNodes is desired number of nodes. If NumNodes is array then
%           number of nodes in calculated graph is not greater than
%           sum(NumNodes). For more details see Several epoch strategies
%           below.
%       ParameterSet is special function to define parameters of
%           algorithm. There are two standard functions:
%               parametersPrincipalCurve.m
%               parametersPrincipalCircle.m
%           Detailed description of each listed function is included into
%           corresponding m files. Description of parameters function
%           development is included into all listed files and to WIKI???
%           If you use parameterFunction with several epoches then see
%           Several epoch strategies below.
%   By default computeElasticPrincipalGraph also prepare tree figures:
%       "Accuracy/Complexity plot" shows the dependence of the normalized
%           geometrical complexity on the number of nodes. This graph is
%           very useful to find optimal number of nodes. You can produce
%           this graph later by function
%           accuracyComplexityPlot(ReportTable) 
%       "PCA view on principal tree" shows the projection of a principal
%           graph onto the first two principal components of the data. In
%           other words this graph is visualisation of data and graph. You
%           can produce this graph later by function PCAView.
%       "MSE and Elastic energy plot" shows the dynamics of data
%           approximation and elastic energy term during principal tree
%           construction. You can produce this graph later by function
%           plotMSDEnergyPlot(ReportTable).
%
%   [NodePositions, Edges, ReportTable] = 
%           computeElasticPrincipalGraph(data, NumNodes, ParameterSet)
%       returns also ReportTable. Report table is table with 16 columns:
%           BARCODE is barcode in form ...S4|S3||N, where N is number of
%               nodes, S3 is number of three stars, S4 (5,...) is number of
%               four (five,...) stars.
%           ENERGY is total energy of graph embedment (ENERGY = MSE + UE +
%               UR)
%           NNODES is number of nodes.
%           NEDGES is number of edges (NEDGES = NNODES - 1).
%           NRIBS is number of two stars (nodes with two otherr connected
%               nodes). 
%           NSTARS is number of stars with 3 and more leaves (nodes
%               connected with central node). 
%           NRAYS2 is sum of rays minus doubled number of nodes.
%           MSE is mean square error or assessment of data approximation
%               quality.
%           MSEP is mean square error of data approximation when the data points are
%               projected to edges rather than graph nodes
%           FVE is fraction of explained variance. This value always
%               between 0 and 1. Greater value means higher quality of
%               data approximation.
%           FVEP is the fraction of explained variance, when the projection
%               of a data point is computed on an edge of the graph and not
%               a node.
%           UE is elastic energy for edges stretching.
%           UR is elastic energy of deviation from harmonicity.
%           URN is UR * nodes 
%           URN2 is UR * nodes^2
%           URSD is standard deviation of UR???
%
%   [...] = computeElasticPrincipalGraph(data,NumNodes,'PARAM1',val1, ...)
%       specifies optional parameter name/value pairs to control the
%       computation. Parameters are: 
%       'Lamdba' is coefficient of elastic stretching of the graph. Coefficient
%           is positive double number. If EP is vector then see Several
%           epoch strategies below.
%       'Mu' is penalty coefficient for deviation from harmonicity.
%           Coefficient is positive double number. If RP is vector then see
%           Several epoch strategies below.
%       'BranchingControl' controls behaviour of the graph for branching
%           The parameter should be a pair of real numbers: [alpha beta]
%           The elasticity of an edge at the graph selection stage 
%           becomes \lambda+\alpha*(k-2), where k
%           is the maximum star order to which the edge belongs. 
%           Thus, higher alpha penalizes appearance of higher order stars.
%           The value of alpha does not affect the optimization of a given 
%           graph structure.
%           The recommended default value for alpha is 0.01.
%           Values of alpha bigger than 0.5 should effectively forbid any
%           branching: thus, forcing construction of a principal curve
%           instead of a tree.
%           The value of beta remains experimental and currently should be
%           set to 1
%       'TrimmingRadius' is robust or trimming radius. To perform non robust
%           algorithm specify TrimmingRadius = 0. 
%       'InitGraph' is structure with two elements: InitNodes and
%           InitEdges. This parameter can be used for continuation of graph
%           construction:
%               % firstly create very harmonic graph
%               [NodePositions, Edges] = ...
%                       computeElasticPrincipalTree(data,10,'RP'=1);
%               % continue for more flexible graph
%               [NodePositions, Edges] = ...
%                       computeElasticPrincipalTree(data,10,'RP'=0.1,...
%                       'InitGraph',struct('InitNodes',NodePositions,...
%                       'InitEdges',Edges));
%       'ReduceDimension' is used for dimensionality reduction by
%           principal components. There are tree possible values:
%               Integer value K. This value must be positive integer being
%                   not greater than m (number of columns in matrix data).
%                   In this case the first K principal components will be
%                   used.
%               Array with two integer values [K1, K2].  These value must
%                   be positive integer being not greater than m (number of
%                   columns in matrix data) and inequality K1 <= K2 has to
%                   be held. In this case principal components from K1 to
%                   K2 inclusive will be used.
%               Array with more than 2 elements. All elements  must be
%                   positive integer being not greater than m (number of 
%                   columns in matrix data). These values are considered as
%                   indices of principal components to use.
%       'Plots' is integer number which specifies set of plots to draw:
%           1 for "Accuracy/Complexity plot" 
%           2 for "PCA view on principal graph"
%           4 for "Metro map layout of the principal tree"
%           8 for "MSE and Elastic energy plot"
%           To specify several plots it is necessary to use sum of listed
%           values. For example, for "Accuracy/Complexity plot" and "Metro
%           map layout of the principal tree" i is necessary to use value
%           5 = 1 (Accuracy/Complexity) + 4 (Metro map).
%
%
%
%Example of use:  
%   load('test_data/iris/iris.mat');
%   [NodePositions, Edges, ReportTable] =...
%       computeElasticPrincipalGraph(table2array(iris(:,2:end)),20);
%
    % Can be removed if all paths are set before.
    setallpaths;

    % Parse optional argumentes
    reduceDimension = 0;
    newDimension = -1;
    drawAccuracyComplexity = true;
    drawPCAView = true;
    drawEnergy = true;
    graphinitialized = 0;
    Lambda = 0.01;
    Mu = 0.1;
    InitStruct = [];
    for i=1:2:length(varargin)
        if strcmpi(varargin{i}, 'ReduceDimension')
            reduceDimension = 1;
            newDimension = varargin{i + 1};
        elseif strcmpi(varargin{i}, 'Plots')
            tmp = uint8(varargin{i + 1});
            drawEnergy = bitand(tmp, 8)>0;
            drawPCAView = bitand(tmp, 2)>0;
            drawAccuracyComplexity = bitand(tmp, 1)>0;
        elseif strcmpi(varargin{i}, 'Lambda')
            Lambda = varargin{i + 1}; 	
        elseif strcmpi(varargin{i}, 'Mu')
            Mu = varargin{i + 1}; 	
        elseif strcmpi(varargin{i}, 'InitGraph')
            InitStruct = varargin{i + 1};
        end
    end

    mv = mean(data);
    data_centered = bsxfun(@minus, data, mv);
    
    if ~isempty(InitStruct)
        % Initialise graph if it is necessary
        np = InitStruct.InitNodes;
        ed = InitStruct.InitEdges;
        if(~isfield(InitStruct,'ElasticMatrix'))
            em = MakeUniformElasticMatrix(ed, Lambda, Mu);
        else
            em = InitStruct.ElasticMatrix;
        end
        np_centered = bsxfun(@minus, np, mv);
        graphinitialized = 1;
    end


    [vglobal, uglobal, explainedVariances] = pca(data_centered);
    if reduceDimension
        % Form index of used PCs
        tmp = length(newDimension);
        if tmp == 1
            indPC = 1:newDimension(1);
        elseif tmp == 2
            indPC = newDimension(1):newDimension(2);
        else
            indPC = newDimension;
        end
        % Calculate variance explained by seleced components.
        perc = sum(explainedVariances(indPC))/sum(explainedVariances)*100;
        display(sprintf('Variance retained in %3.0f dimensions: %2.2f%%',...
            (length(indPC)),perc));
        data_centered = uglobal(:,indPC);
        if graphinitialized
            np_centered = np_centered*vglobal(:,indPC);
        end
    else
        indPC = 1:size(data,2);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Computing the graph

    if ~graphinitialized
        [NodePositions, ElasticMatrix, ReportTable] =...
            ElPrincGraph(data_centered, NumNodes, Lambda, Mu, varargin{:});
    else
        [NodePositions, ElasticMatrix, ReportTable] =...
            ElPrincGraph(data_centered, NumNodes, Lambda, Mu,...
            'InitNodePositions', np_centered, 'InitElasticMatrix', em, varargin{:});
    end
    
    [row, col] = find(triu(ElasticMatrix, 1));
    Edges = [row, col];

    %%%%%%%%%%%%%%%%%%%%%%%% Preparing the output arguments
    if reduceDimension
        %Project nodes back into the initial, non-reduced space
        NodePositions = NodePositions * vglobal(:, indPC)';
    end
    NodePositions = bsxfun(@plus, NodePositions, mv);

    %%%%%%%%%%%%%%%%%%%%%%%%  Plots of MSE, elastic energy optimization
    if drawEnergy
        plotMSDEnergyPlot(ReportTable,explainedVariances);
        set(gcf, 'Position', [5   102   499   171]);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%  Accuracy/Complexity plot
    if drawAccuracyComplexity
        accuracyComplexityPlot(ReportTable);
        set(gcf, 'Position', [8   361   499   275]);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%% Show principal component view on principal
    %%%%%%%%%%%%%%%%%%%%%%%% tree and the data
    % This method use two the first selected component in the list of used
    % components (indPC)
    
    if drawPCAView
        PCAView( NodePositions, Edges, data,...
            vglobal(:, indPC(1)), vglobal(:, indPC(2)),... 
            explainedVariances(indPC(1)) / sum(explainedVariances),...
            explainedVariances(indPC(2)) / sum(explainedVariances),...
            varargin{:});
        set(gcf, 'Position', [511   156   510   413]);
    end
    
    drawnow;
end