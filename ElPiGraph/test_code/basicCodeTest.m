%%%%%%%%%%%%%%%%% Test 1
display('==================== Test 1 ========================');

display('Constructing principal tree for iris data');

data = load('test_data/iris/iris.data');
[NodePositions,Edges,ReportTable] = computeElasticPrincipalGraph(data,20);

pause(2);
close all;

display('Constructing principal curve for circle example');

cc = load('test_data/circle/simple_circle.data');
computeElasticPrincipalCurve(cc,30);

pause(2);
close all;

display('Constructing principal circle for circle example');

cc = load('test_data/circle/simple_circle.data');
computeElasticPrincipalCircle(cc,30);

pause(2);
close all;


display('Changing elasticity moduli (more rigid graph)');

computeElasticPrincipalGraph(data,20,'Lambda',0.1,'Mu',0.5);

pause(2);
close all;

display('Robust version, trimming radius set to 0.7');

computeElasticPrincipalGraph(data,40,'TrimmingRadius',0.7);

pause(2);
close all;

display('Constructing principal tree in the first two PC dimension');

computeElasticPrincipalGraph(data,30,'ReduceDimension',2);


pause(2);
close all;

display('Constructing principal tree after subtracting the first PC');

computeElasticPrincipalGraph(data,30,'ReduceDimension',[2 size(data,2)]);

pause(2);
close all;

display('Constructing principal tree in two epochs, rigid and soft');

[NodePositions,Edges,ReportTable] = computeElasticPrincipalGraph(data,20,'Lambda',0.1,'Mu',1);
close all;     
computeElasticPrincipalGraph(data,30,'InitGraph',struct('InitNodes',NodePositions,'InitEdges',Edges),'Lambda',0.001,'Mu',0.01);


pause(2);
close all;

display('Showing only PCA view plot');

[NodePositions,Edges,ReportTable] = computeElasticPrincipalGraph(data,20,'Plots',2);

pause(2);
close all;

display('Comparing global and local optimization on a relatively large graph and dataset');

data = load('test_data/tree23/tree23_inflated.data');
figure;
display('Global optimization:');
tic; [NodePositions,Edges,ReportTable] = computeElasticPrincipalGraph(data,50,'Plots',2,'verbose',0); toc;
display('Local optimization:');
figure;
tic; [NodePositions,Edges,ReportTable] = computeElasticPrincipalGraph(data,50,'Plots',2,'LocalSearch',2,'verbose',0); toc;

pause(2);
close all;

display('Constructing principal tree without shrinking grammar');

data = load('test_data/tree23/tree23_inflated.data');
tic;
[NodePositions,Edges,ReportTable] = computeElasticPrincipalGraph(data,50,'GrowGrammars',[{'bisectedge','addnode2node'}],'ShrinkGrammars',[]);
toc;

%%%%%%%%%%%%%%%%% Test 2
display('==================== Test 2 ========================');

display('Simple test for a branching dataset');

data = load('test_data\tree23\tree23.data'); 
plot(data(:,1),data(:,2),'ko','MarkerSize',5);
[NodePositions,Edges] = computeElasticPrincipalGraph(data,30);

pause(2);
figure;

display('Creating Metro Map layout');

 % compute the number of data points projected into each node
 partition = PartitionData(data,NodePositions,100000,sum(data.^2,2));
 NodeSizes = histc(partition,[min(partition):max(partition)]);
 % draw the metro map layout, visualize the number of data points projected
 drawMetroMap(NodePositions,Edges,'NodeSizes',NodeSizes)

pause(2);
close all;

display('Simple test for a branching dataset with noise (~20%% of noise points)');

data_noise = load('test_data\tree23\tree23_noise.data'); 
plot(data_noise(:,1),data_noise(:,2),'ko','MarkerSize',5);

pause(1);

[NodePositions,Edges] = computeElasticPrincipalGraph(data_noise,30);

pause(2);
close all;

display('Producing robust graph');

[NodePositions,Edges] = computeRobustElasticPrincipalGraph(data_noise,40,0.4,'Lambda',0.02);

pause(2);
close all;

display('Adding 200%% of noisy data points');

amount_of_noise = 200;
data_noise = add_uniform_background_noise(data,amount_of_noise);
[NodePositions,Edges] = computeRobustElasticPrincipalGraph(data_noise,30,0.3,'Lambda',0.02,'Plots',0);
plot(data_noise(:,1),data_noise(:,2),'k.','MarkerSize',5); hold on; 
plot(data(:,1),data(:,2),'r.','MarkerSize',5); 
drawGraph2D(NodePositions,Edges); 
title(sprintf('Percentage of noise points = %i%%',amount_of_noise));

pause(2);
close all;


%%%%%%%%%%%%%%%%% Test 3
display('==================== Test 3 ========================');

display('Creating elastic principal tree');

iristable = importdata('test_data\iris\iris.txt'); 
iris = iristable.data;
irisz = zscore(iris);
pointlabels = (iristable.textdata(2:end,1))';

[NodePositions,Edges] = computeElasticPrincipalGraph(irisz,20);


pause(2);
close all;

display('Showing Metro Map Layout');


NodesMM = computeMetroMapLayout(NodePositions,Edges);
close all; 
drawGraph2D(NodesMM,Edges,'ShowClusterNumbers',0); 
lcm = drawPieChartsMetroMap(NodePositions,Edges,irisz,pointlabels,NodesMM);


pause(2);
close all;

display('Showing Metro Map Layout, specifying colors for pie-charts');


lcm('Iris-setosa') = [1 1 0]; % Setosa will be yellow
lcm('Iris-versicolor') = [0 1 1]; % Versicolor will be light blue
lcm('Iris-virginica') = [1 0 1]; % Virginica will be purple
close all; 
drawGraph2D(NodesMM,Edges,'ShowClusterNumbers',0); 
drawPieChartsMetroMap(NodePositions,Edges,irisz,pointlabels,NodesMM,'LabelColorMap',lcm);

pause(2);
close all;
