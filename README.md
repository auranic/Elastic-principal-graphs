# Elastic principal graphs
Matlab implementation of Elastic Principal Graphs (ElPiGraph) method. 

For the basic self-contained formal description of ElPiGraph, read https://github.com/auranic/Elastic-principal-graphs/blob/master/ElPiGraph_Methods.pdf. 

Principal graphs are graphs that are embedded into a high-dimensional space and minimize the distance to the data points, while maximizing some regular properties.

Elastic principal graphs are based on minimization of the energy potential containing three parts :

### U = MSE + \lambda UE + \mu UR

where MSE is the mean squared error of data approximation, UR - is the sum of squared edge lengths, UR is a term minimizing the deviation of the principal graph embedment from harmonicity (generalization of linearity), \lambda ,\mu  are coefficients of regularization.

The structure of the graph is computed by an optimal application of a sequence of graph transformations, using operations from predefined graph grammar.
The simplest graph grammar "Bisect an edge", "Add a node to a node" leads to construction of a principal tree.

Read [wiki](https://github.com/auranic/Elastic-principal-graphs/wiki) of this repository for more detailed information on the algorithm
and examples of its application.

Go to [R implementation of elastic principal graphs interface](https://github.com/Albluca/rpgraph) created by [Luca Albergante](https://github.com/Albluca).

For more details of elastic principal graph theory read:

[Gorban AN and Zinovyev AY. Principal Graphs and Manifolds. In Handbook of Research on Machine Learning Applications and Trends: Algorithms, Methods and Techniques (eds. Olivas E.S., Guererro J.D.M., Sober M.M., Benedito J.R.M., Lopes A.J.S.). Information Science Reference, September 4, 2009.](https://arxiv.org/ftp/arxiv/papers/0809/0809.0490.pdf)

[Gorban A., Sumner N., Zinovyev A. Topological grammars for data approximation. 2007. Applied Mathematics Letters 20(4), 382-386.](http://arxiv.org/pdf/cs/0603090v2.pdf)

[Gorban A., Mirkes E., Zinovyev A. Robust principal graphs for data approximation. Archives of Data Science 2(1):1:16, 2017.](href="http://arxiv.org/abs/1603.06828)

[Gorban A.N., Zinovyev A. 2010. Principal manifolds and graphs in practice: from molecular biology to dynamical systems. Int J Neural Syst 20(3):219-32.](href="http://arxiv.org/pdf/1001.1122.pdf)

# Organization of the code

## Folders:

core_algorithm 		- contains the core MATLAB code of the algorithm (self-contained)

core_algorithm_java	- old MATLAB wrapper of the Java code for ElPiGraph

docs			- some documentation on the method and the code

examples		- some example applications of the method

simulations		- code generating synthetic datasets (e.g., with known branching topology)

test_code		- testing critical parts of the code (not needed for the package use)

test_data		- example datasets 

utils			- utility functions for manipulating data and the graph (e.g., projection function)

visualization		- utility functions used for visualizing the results of the method application (such as applying metro map layout of the principal tree)

## Functions in the root folder:

computeElasticPrincipalCircle.m		- computes closed elastic principal curve with a given number of nodes

computeElasticPrincipalCurve.m          - computes elastic principal curve with a given number of nodes

computeElasticPrincipalGraph.m		- computes elastic principal graph for a dataset with a given number of nodes and defined set of grammars (principal tree grammar by default)

setallpaths.m				- set all paths to other folders (needed if some functions are called directly, bypassing the root folder "computeElasticPrincipalXXX.m" functions)

 
