# Elastic principal graphs
Matlab implementation of the interface to the  Elastic Principal Graphs method

Principal graphs are graphs that are embedded into a high-dimensional space and minimize the distance to the data points, while maximizing some regular properties.

Elastic principal graphs are based on minimization of the energy potential containing three parts :

### U = MSE + ep UE + rp UR

where MSE is the mean squared error of data approximation, UR - is the sum of squared edge lengths, UR is a term minimizing the deviation of the principal graph embedment from harmonicity (generalization of linearity), ep,rp are coefficients of regularization.

The structure of the graph is computed by an optimal application of a sequence of graph transformations, using operations from predefined graph grammar.
The simplest graph grammar "Bisect an edge", "Add a node to a node" leads to construction of a principal tree.

Read wiki of this repository for much more detailed information on the algorithm.

For more details of elastic principal graph theory consult the bibliography:

1) Gorban A.N., Zinovyev A. 2010. Principal manifolds and graphs in practice: from molecular biology to dynamical systems. Int J Neural Syst 20(3):219-32.

2) Gorban AN and Zinovyev AY. Principal Graphs and Manifolds. In Handbook of Research on Machine Learning Applications and Trends: Algorithms, Methods and Techniques (eds. Olivas E.S., Guererro J.D.M., Sober M.M., Benedito J.R.M., Lopes A.J.S.). Information Science Reference, September 4, 2009.

3) Gorban A, Sumner N, Zinovyev A. Beyond The Concept of Manifolds: Principal Trees, Metro Maps, and Elastic Cubic Complexes. 2008. Lecture Notes in Computational Science and Engeneering 58: 223-240.

4) Gorban AN, Sumner NR, Zinovyev AY. Topological grammars for data approximation. Applied Mathematics Letters 20 (4), 382-386.

5) Gorban A., Zinovyev A. Elastic Principal Graphs and Manifolds and their Practical Applications. 2005. Computing 75,359 -379
