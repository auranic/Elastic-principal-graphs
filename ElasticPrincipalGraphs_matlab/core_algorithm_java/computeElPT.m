function [NodePositions, Edges, ReportTable, cpg] = ...
                                    computeElPT(data, NumNodes, varargin)
% computeElPT is service function to calculate principal graph. It is
% preferable to use functions computeElasticPrincipalTree,
% computeElasticPrincipalGraph or makeAnimation.
%
%Usage
%   [NodePositions, Edges] = computeElPT(data,NumNodes)
%   returns 
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
%
%   [NodePositions, Edges, ReportTable] = computeElPT(data,NumNodes)
%       returns also ReportTable. Report table is table with 17 columns:
%           STEP is number of iteration
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
%           MSEP is project to edge ???
%           FVE is fraction of explained variance. This value always
%               between 0 and 1. Greater value means higher quality of
%               data approximation.
%           FVEP is ???
%           UE is elastic energy for edges stretching.
%           UR is elastic energy of deviation from harmonicity.
%           URN is UR * nodes 
%           URN2 is UR * nodes^2
%           URSD is standard deviation of UR???
%
%   [NodePositions, Edges, ReportTable, cpg] = computeElPT(data,NumNodes)
%       returns also vdaoengine.analysis.grammars.ComputePrincipalGraph
%           Java container object cpg with various service functions.
%           Description of this functions can be found in ???
%
%   [...] = computeElasticPrincipalTree(data,NumNodes,'PARAM1',val1, ...)
%       specifies optional parameter name/value pairs to control the
%       computation. Parameters are: 
%       'EP' is coefficient of elastic stretching of the graph. Coefficient
%           is positive double number. If EP is vector then see Several
%           epoch strategies below.
%       'RP' is penalty coefficient for deviation from harmonicity.
%           Coefficient is positive double number. If RP is vector then see
%           Several epoch strategies below.
%       'ParameterSet' is special function to define parameters of
%           algorithm. There are tree standard functions:
%               parametersDefaultPrincipalTree.m
%               parametersRobustElasticPrincipalTree.m
%               parametersTwoStageRobustElasticPrincipalTree.m
%           Default value is parametersDefaultPrincipalTree. Detailed
%           description of each listed function is included into
%           corresponding m files. Description of parameters function
%           development is included into all listed files and to WIKI???
%           If you use parameterFunction with several epoches (for example,
%           parametersTwoStageRobustElasticPrincipalTree.m) then see Several
%           epoch strategies below.
%       'TrimRadius' is robust or trimming radius. To perform non robust
%           algorithm specify TrimRadius = 0. If TrimRadius is array then
%           each value is used for separate epoch. see Several epoch
%           strategies below. 
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
%
%   Several epoch strategies. This software call epoch fragment of
%       algorithm which is processed with the same set of patrametres.
%       There are several ways to create strategy with several epochs. The
%       simplest is use parameter function
%       parametersTwoStageRobustElasticPrincipalTree, which provide two
%       epochs: the first epoch is performed as non-robust and the second
%       epoch as robust. Sometimes it is very useful initially create
%       graph with high EP and RP and then decrease this parameters. To do
%       this without special development of parameter function you can use
%       four parameters to customise algorithm: 
%           NumNodes is desired number of nodes.
%           'EP' is coefficient of elastic stretching of the graph. 
%           'RP' is penalty coefficient for deviation from harmonicity.
%           'TrimRadius' is robust or trimming radius.
%       Final strategy will be formed by algorithm:
%       Number of epoch in strategy is  
%           nEpochs = max([length(RP), length(EP), length(NumNodes),...
%                   length(TrimRadius), length(parameters.epochs)]);
%           for k = 1:nEpochs
%               %Create new epoch object.
%               % if we have corresponding element in parameters.epochs
%               % then we use it otherwise use the last epoch
%               if k <= length(parameters.epochs)
%                   copy parameters.epochs(k);
%               else
%                   copy parameters.epochs(end);
%               end
%               % Check modification by other parameters
%               % if we have corresponding element in EP then use it
%               if k <= length(EP)
%                   if EP(k)>0
%                       epoch.EP = EP(k);
%                   end
%               end
%               % if we have corresponding element in RP then use it
%               if k <= length(RP)
%                   if RP(k)>0
%                       epoch.RP = RP(k);
%                   end
%               end
%               % if we have corresponding element in NumNodes then use it
%               % oterwise use the last element
%               if k <= length(RP)
%                   epoch.numberOfIterations = NumNodes(k);
%               else
%                   epoch.numberOfIterations = NumNodes(length(NumNodes));
%               end
%               % if we have corresponding element in TrimRadius then use it
%               if k <= length(TrimRadius)
%                   if TrimRadius(k)>0
%                       epoch.robust = true;
%                   end
%                   epoch.trimradius = TrimRadius(k);
%               end
%               % Add new epoch to 
%               config.epochs.add(epoch);
%           end
%
%Examples can be found in the following functions:
%       computeElasticPrincipalTree,
%       computeElasticPrincipalGraph,
%       makeAnimation.
%
                                
    % Elasticity module stretching
    EP = -1;
    % Elasticity module bending
    RP = -1;
    % Trimming radius
    TrimRadius = -1;
    % Graph initialisation
    initGraph = 0;
    % Function to handle parameters
    parameterfunction_handle = @parametersDefaultPrincipalTree;
    % Parse input parameters
    for i=1:2:length(varargin)
        if strcmpi(varargin{i},'EP')
            EP = varargin{i+1};
        elseif strcmpi(varargin{i},'RP')
            RP = varargin{i+1};
        elseif strcmpi(varargin{i},'ParameterSet')
            parameterfunction_handle = varargin{i+1};
        elseif strcmpi(varargin{i},'TrimRadius')
            TrimRadius = varargin{i+1};
        elseif strcmpi(varargin{i},'InitGraph')
            initGraph = 1;
            tmp = varargin{i+1};
            InitNodePositions = tmp.InitNodes
            InitEdges = tmp.InitEdges;
        end
    end

    % interface to java calculator
    javaclasspath({'VDAOEngine.jar'});
    javaclasspath({'core_algorithm_java/VDAOEngine.jar'});

    % Get parameters from parameter function;
    parameters = parameterfunction_handle();
    % Form configuration object
    config = vdaoengine.analysis.grammars.ConfigFile;
    % Epecify parameters of configuration
    config.algtype = parameters.algtype;
    config.stretchInitCoeffs(1) = parameters.initstretch1;
    config.initStrategy = parameters.initalgorithm;
    % Work with several epochs in configuration
    % Define number of epochs
    nEpochs = max([length(RP),length(EP), length(NumNodes),...
        length(TrimRadius), length(parameters.epochs)]);
    for k = 1:nEpochs
        % Create new java epoch and fill it.
        epoch = vdaoengine.analysis.elmap.ElmapAlgorithmEpoch;
        % if we have corresponding element in parameters.epochs then we use
        % it otherwise use the last epoch
        if k <= length(parameters.epochs)
            i = k;
        else
            i = length(parameters.epochs);
        end
        epoch.grammarType = parameters.epochs(i).grammartype;
        epoch.EP = parameters.epochs(i).ep;
        epoch.RP = parameters.epochs(i).rp;
        epoch.maxNumberOfIterationsForSLAU = parameters.epochs(i).numiterSLAU;
        epoch.epsConvergence = parameters.epochs(i).eps;
        epoch.epsConvergenceSLAU = parameters.epochs(i).epsSLAU;
        if isfield(parameters.epochs(i),'robust')
            epoch.robust = parameters.epochs(i).robust;
        end
        if isfield(parameters.epochs(i),'trimradius')
            epoch.trimradius = parameters.epochs(i).trimradius;
        end
        % Check modification by other parameters
        % if we have corresponding element in EP then use it
        if k <= length(EP)
            if EP(k)>0
                epoch.EP = EP(k);
            end
        end
        % if we have corresponding element in RP then use it
        if k <= length(RP)
            if RP(k)>0
                epoch.RP = RP(k);
            end
        end
        % if we have corresponding element in NumNodes then use it oterwise
        % use the last element
        if k <= length(RP)
            epoch.numberOfIterations = NumNodes(k);
        else
            epoch.numberOfIterations = NumNodes(length(NumNodes));
        end
        % if we have corresponding element in TrimRadius then use it
        if k <= length(TrimRadius)
            if TrimRadius(k)>0
                epoch.robust = true;
            end
            epoch.trimradius = TrimRadius(k);
        end
        % Add new epoch to 
        config.epochs.add(epoch);
    end
    % Create object for calculation and set data and configuretion
    cpg = vdaoengine.analysis.grammars.ComputePrincipalGraph;
    cpg.config = config;
    cpg.setDataSetAsMassif(data);
    
    if initGraph > 0
        cpg.config.initStrategy = -1;
        cpg.setPrimitiveGraphByNodesAndEdges(InitNodePositions,InitEdges);
    end
    % Calculations
    report = cpg.compute();
    % Convert report to usable format
    fn = tempname;
    fid = fopen(fn,'w');
    fprintf(fid, '%s', char(report));
    fclose(fid);
    ReportTable = readtable(fn,'Delimiter','\t');
    NodePositions = cpg.graph.getNodePositions();
    Edges = cpg.graph.getEdgeTable();
end