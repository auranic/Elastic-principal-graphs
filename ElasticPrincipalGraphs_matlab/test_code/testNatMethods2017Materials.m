mat = load('./test_data/DDRTree2017Materials/tree_300.mat');

X = mat.X; 

Npoints = size(X,1)
dim = size(X,2)


maxNumNodes = 50;

[np,ed,ReportTable] = computeElasticPrincipalGraph(X,maxNumNodes,'Lambda',0.03); 

%[np,ed,ReportTable] = computeElasticPrincipalGraph(X,maxNumNodes,'Lambda',0.01); 
