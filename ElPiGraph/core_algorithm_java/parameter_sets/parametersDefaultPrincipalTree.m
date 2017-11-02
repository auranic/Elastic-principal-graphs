function params = parametersDefaultPrincipalTree()
%parametersDefaultPrincipalTree prepare set of parameters for Elastic
%Principal Grapph package.
%
%Parameters function does not have input arguments.
%Parameter function returns structure params with following set of fields
%   algtype (required) is type of algorithm. For this package there is one
%       possible value only:
%           'grammar' is algorithm based on the graph grammars.
%   initstretch1 (required) is initial coefficient of stretching in created
%       graph. This parameter is positive real number. Default value is 1.
%   initalgorithm (required) is number of algorithm to ce]reate initial
%       graph. List of acceptable values is:
%           0 - initialization by a piece of first principal component.
%           1 - initialization by a piece of first principal component with
%               alignment to trimming radius. ???
%           2 - initialization by a closed reactangle spanned by first two
%               PCs.
%   epochs (required) is array of epochs of graph fitting. At least one
%       element must be specified. User can specify usage of several
%       epoches without redevelopment of parameter function (see Several
%       epoch strategies in computeElPT, computeElasticPrincipalTree, or
%       computeElasticPrincipalGraphEach). If user specified more epochs
%       than parameter function then the last epoch from parameter function
%       will be used for remaining epoches. Each epoch is structure with
%       following fields:  
%       grammartype (required) is type of grammar to use in this epoch.
%           Grammar type can have following values:
%           'treeWithPruning' is grammar which use 2 steps of tree growing
%               with rules "Add a node to a node" and "Bisect an edge"
%               and one step of tree pruning with rules "Remove terminal
%               node" and "Glue two connected nodes".
%           'circle' is is the simplest grammar and contains rule "Bisect
%               an edge" only.
%           'curve' is grammar with two riles "Add a node to a terminal
%               node" and "Bisect an edge".
%       ep (required) is coefficient of elastic stretching of the graph.
%           Coefficient is positive double number. This value can be
%           redefined by parameter 'EP' in functions which uses parameter
%           function. (see Several epoch strategies in computeElPT,
%           computeElasticPrincipalTree, or computeElasticPrincipalGraphEach).
%       rp (required) is penalty coefficient for deviation from harmonicity.
%           Coefficient is positive double number. This value can be
%           redefined by parameter 'EP' in functions which uses parameter
%           function. (see Several epoch strategies in computeElPT,
%           computeElasticPrincipalTree, or computeElasticPrincipalGraphEach).
%       numiterSLAU (required) is number of iteration of Gauss–Seidel
%           method to solve system of linear algebraic equations.
%           Recommended value is 1000.
%       epsSLAU (required) is required accuracy of convergence of
%           Gauss–Seidel method. Recommended value is 0.01.
%       eps (required) is required accuracy of convergence of node position
%           optimisation. Recommended value is 0.01. ??? Is my guess correct?
%       robust (optional) is true for robust graph and false otherwise.
%           Default value is false. This value can be redefined by
%           parameter 'TrimRadius' in functions which uses parameter 
%           function. (see Several epoch strategies in computeElPT,
%           computeElasticPrincipalTree, or computeElasticPrincipalGraphEach).
%       trimradius (optional) is trimming radius for robust graph.
%           Default value is 0.1 (??? 0.1 of what???). This value can be redefined by
%           parameter 'TrimRadius' in functions which uses parameter 
%           function. (see Several epoch strategies in computeElPT,
%           computeElasticPrincipalTree, or computeElasticPrincipalGraphEach).
%       
%There are several depricated fields of epoch structure which ARE NOT USED
%in current version of package but are described in examples of parameter
%functions : 
%   id, minimize, numiter.
%

    params.algtype = 'grammar';

    params.initstretch1 = 1;

    % initialization by a piece of first principal component
    params.initalgorithm = 0;

    % universal parameters (should be rarely modified)
    params.epochs(1).id = 1;
    params.epochs(1).minimize = 1;
    params.epochs(1).eps = 0.01;
    params.epochs(1).epsSLAU = 0.01;
    params.epochs(1).numiterSLAU = 1000;

    % what grammar to use
    params.epochs(1).grammartype = 'treeWithPruning';

    % elasticity parameters
    params.epochs(1).ep = 0.01;
    params.epochs(1).rp = 0.1;

    % default number of nodes
    params.epochs(1).numiter = 20;
end