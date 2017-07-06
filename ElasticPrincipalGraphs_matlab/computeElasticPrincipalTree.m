function [NodePositions, Edges, ReportTable, cpg] =...
    computeElasticPrincipalTree(data, NumNodes, varargin)
%computeElasticPrincipalTree calculate elastic principal tree which is one
%of the simplest and very useful type of principal graphs.
%Theory and some examples of computeElasticPrincipalTree usage can be found
%in https://github.com/auranic/Elastic-principal-graphs/wiki/Basic-use-of-Elastic-Principal-Graphs-Matlab-package
%
%Usage
%   [NodePositions, Edges] = computeElasticPrincipalTree(data,NumNodes)
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
%   By default computeElasticPrincipalTree also prepare four figures:
%       "Accuracy/Complexity plot" shows the dependence of the normalized
%           geometrical complexity on the number of nodes. This graph is
%           very useful to find optimal number of nodes. You can produce
%           this graph later by function
%           accuracyComplexityPlot(ReportTable) 
%       "PCA view on principal tree" shows the projection of a principal
%           graph onto the first two principal components of the data. In
%           other words this graph is visualisation of data and graph. You
%           can produce this graph later by function PCAView.
%       "Metro map layout of the principal tree" presents 'topological
%           skeleton' of graph. In the PCA view of principal tree some
%           nodes can be hidden by other nodes but in metro map all nodes
%           always visible. You can produce this graph later by function
%           drawMetroMap(NodePositions,Edges, 'NodeSizes',...
%               cpg.graph.countNumberOfPointsProjected(cpg.dataset) + 1); 
%       "MSE and Elastic energy plot" shows the dynamics of data
%           approximation and elastic energy term during principal tree
%           construction. You can produce this graph later by function
%           plotMSDEnergyPlot(ReportTable).
%
%   [NodePositions, Edges, ReportTable] =
%                              computeElasticPrincipalTree(data,NumNodes)
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
%   [NodePositions, Edges, ReportTable, cpg] = 
%                               computeElasticPrincipalTree(data,NumNodes)
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
%       'Reduce dimension' is indicator of dimensionality reduction by
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
%       'Graphs' is integer number which specifies set of graphs to draw:
%           1 for "Accuracy/Complexity plot" 
%           2 for "PCA view on principal tree"
%           4 for "Metro map layout of the principal tree"
%           8 for "MSE and Elastic energy plot"
%           To specify several graphs it is necessary to use sum of listed
%           values. For example, for "Accuracy/Complexity plot" and "Metro
%           map layout of the principal tree" i is necessary to use value
%           5 = 1 (Accuracy/Complexity) + 4 (Metro map).
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
%
%
%Example:
%   load('test_data/iris/iris.mat');
%   [NodePositions, Edges, ReportTable, cpg, mml] =...
%       computeElasticPrincipalTree(table2array(iris(:,2:end)),20);
%

    % Parse optional argumentes
    reduceDimension = 0;
    newDimension = -1;
    drawAccuracyComplexity = true;
    drawPCAView = true;
    drawMetroMaps = true;
    drawEnergy = true;

    for i=1:2:length(varargin)
        if strcmpi(varargin{i},'ReduceDimension')
            reduceDimension = 1;
            newDimension = varargin{i+1};
        elseif strcmpi(varargin{i},'Graphs')
            tmp = uint8(varargin{i+1});
            drawEnergy = bitand(tmp,8)>0;
            drawMetroMaps = bitand(tmp,4)>0;
            drawPCAView = bitand(tmp,2)>0;
            drawAccuracyComplexity = bitand(tmp,1)>0;
        end
    end

    mv = mean(data);
    data_centered = bsxfun(@minus,data,mv);

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
    else
        indPC = 1:size(data,2);
    end

	% Calculate principal tree
    [NodePositions, Edges, ReportTable, cpg] = ...
        computeElPT(data_centered,NumNodes,varargin{:});

    %%%%%%%%%%%%%%%%%%%%%%%% Preparing the output arguments
    
    if reduceDimension
        %Project nodes back into the initial, non-reduced space
        NodePositions = NodePositions*vglobal(:,indPC)';
    end
    NodePositions = bsxfun(@plus,NodePositions,mv);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%  Plots of MSE, elastic energy optimization
    if drawEnergy
        plotMSDEnergyPlot(ReportTable,explainedVariances);
        set(gcf,'Position',[11   338   612   300]);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%  Accuracy/Complexity plot
    if drawAccuracyComplexity
        accuracyComplexityPlot(ReportTable);
        set(gcf,'Position',[11   665   612   316]);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%% Show principal component view on principal
    %%%%%%%%%%%%%%%%%%%%%%%% tree and the data
    % This method use two the first selected component in the list of used
    % components (indPC)
    
    if drawPCAView
        PCAView( NodePositions, Edges, data, cpg,...
            vglobal(:,indPC(1)), vglobal(:,indPC(2)),... 
            explainedVariances(indPC(1))/sum(explainedVariances),...
            explainedVariances(indPC(2))/sum(explainedVariances));
        set(gcf,'Position',[11   665   612   316]);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%% Producing metro map layout
    if drawMetroMaps
        drawMetroMap(NodePositions,Edges, 'NodeSizes',...
            cpg.graph.countNumberOfPointsProjected(cpg.dataset) + 1);
        set(gcf,'Position',[841   239   695   643]);
    end
end