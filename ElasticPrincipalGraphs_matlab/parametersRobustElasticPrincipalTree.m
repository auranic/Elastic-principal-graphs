function params = parametersTwoStageRobustElasticPrincipalTree()

params.algtype = 'grammar';

params.initstretch1 = 1;

% initialization by a piece of first principal component
params.initalgorithm = 1;

% universal parameters (should be rarely modified)
params.epochs(1).id = 1;
params.epochs(1).minimize = 1;
params.epochs(1).eps = 0.01;
params.epochs(1).epsSLAU = 0.01;
params.epochs(1).numiterSLAU = 1000;

% what grammar to use
params.epochs(1).grammartype = 'treeWithPruning';

% elasticity parameters
params.epochs(1).ep = 0.001;
params.epochs(1).rp = 0.1;

% default number of nodes
params.epochs(1).numiter = 30;

% robustness parameters
params.epochs(1).robust = true;
params.epochs(1).trimradius = 0.1;