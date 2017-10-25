setallpaths;


global numberOfLocalPoints;
global numberOfFits;
global numberOfFullFits;
global numberOfGraphNodes;
global FractionOfGraphNodes;
global TimeForFitting;
global CountNumberOfIterations;


numberOfLocalPoints=[];
numberOfFits=0;
numberOfFullFits=0;
numberOfGraphNodes=0;
FractionOfGraphNodes=[];
TimeForFitting=[];
NumberOfIterations=[];

%X = load('./test_data/iris/iris.data');
X = load('./test_data/tree23/tree23.data');
%X = load('./test_data/circle/simple_circle.data');
%X = load('./test_data/mosaic/mosaic.txt');
%X = load('C:/Datas/ElGraph_Matlab/competitors/MERLOT/cells.txt');
%X = load('C:/Datas/ElGraph_Matlab/competitors/PAMI2016/toy/tree_600.txt');

%X = load('./test_data/development/guo2010.data');

%X = zscore(X);

% inflate the number of points
inflationFactor = 100;
X1 = zeros(size(X,1)*inflationFactor,size(X,2));
k=1;
STDV = std(X);
if 1
for i=1:size(X,1)
    for j=1:inflationFactor
        r = rand(1,size(X,2));
        p = 0.5*STDV.*r;
        X1(k,:) = X(i,:)+p;
        k=k+1;
    end
end
X = X1;
end

maxNumNodes = 20;

Npoints = size(X,1)
dim = size(X,2)

%tic; computeElasticPrincipalTree(X,maxNumNodes,'EP',0.1); toc;


   profile on
   %tic; [np,em] = ElPrincTree_matlab(X,maxNumNodes,0.01,0.1); toc;
   %tic; [np,ed,ReportTable] = computeElasticPrincipalGraph_java(X,maxNumNodes,@parametersPrincipalCircle,'RP',0.0001); toc;
   close all;
   tic; [np,ed,ReportTable] = computeElasticPrincipalGraph(X,maxNumNodes); toc;
   %tic; [np,ed,ReportTable] = computeElasticPrincipalGraph(X,maxNumNodes,'LocalSearch',2); toc;
   profile viewer
   %tic; [np,ed,ReportTable] = computeElasticPrincipalCurve(X,maxNumNodes,'Mu',10); toc;
   %tic; [np,ed,ReportTable] = computeElasticPrincipalCircle(X,maxNumNodes,'Mu',0.0001); toc;
   %tic; [np,ed,ReportTable] = computeElasticPrincipalGraph(X,maxNumNodes,'GrowGrammars',[{'bisectedge'}],'ShrinkGrammars',[]); toc;
   em = MakeUniformElasticMatrix(ed,0.01,0.1);
   %figure; plot(X(:,1),X(:,2),'b.'); hold on; drawGraph2D(np,ed);
   %[v,u,s] = pca(X);
   %figure; plot(u(:,1),u(:,2),'b.'); hold on; drawGraph2D(np*v,ed);
    
   
   display(sprintf('Average number of data points in the fit = %2.2f',mean(numberOfLocalPoints)));
   display(sprintf('Average number of graph nodes in the fit = %2.2f',numberOfGraphNodes/numberOfFits));
   display(sprintf('Average fraction of graph nodes = %2.2f',mean(FractionOfGraphNodes)));
   display(sprintf('Average number of iterations = %2.2f',mean(CountNumberOfIterations)));
   display(sprintf('Average time for fitting = %2.2e',mean(TimeForFitting)));
   display(sprintf('Number of fits = %i',numberOfFits));
   display(sprintf('Number of full fits = %i',numberOfFullFits));
   
   figure;
   hist(FractionOfGraphNodes,50); title('FractionOfGraphNodes');
   figure;
   hist(numberOfLocalPoints,50); title('numberOfLocalPoints');
   figure;
   hist(TimeForFitting,50); title('TimeForFitting');
  
   
  