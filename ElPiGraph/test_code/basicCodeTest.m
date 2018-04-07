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

display('Constructing principal graph for an example with a loop and branch');

cc = load('test_data/circle/loop_branch.data');
[np0,ed0] = computeElasticPrincipalCircle(cc,6,'Plots',0);
[np,ed] = computeElasticPrincipalGraph(cc,30,'InitGraph',struct('InitNodes',np0,'InitEdges',ed0),'BranchingControls',[0.01,1],'Plots',2);

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
close all; drawnow;

display('Comparing global and local optimization on a relatively large graph and dataset');

data = load('test_data/tree23/tree23_inflated.data');

display(sprintf('Number of data points = %i',size(data,1)));

display('Global optimization:');
tic; [NodePositions,Edges,ReportTable] = computeElasticPrincipalGraph(data,50,'Plots',2,'verbose',0); toc;
ln = ones(1,max(max(Edges)));
le = ones(1,size(Edges,1));
display(sprintf('Barcode for the tree = %s',getPrimitiveGraphStructureBarCode(Encode2ElasticMatrix(Edges,le,ln))));
display('Local optimization:');
%figure;
tic; [NodePositions,Edges,ReportTable] = computeElasticPrincipalGraph(data,50,'Plots',2,'LocalSearch',2,'verbose',0); toc;
ln = ones(1,max(max(Edges)));
le = ones(1,size(Edges,1));
display(sprintf('Barcode for the tree = %s',getPrimitiveGraphStructureBarCode(Encode2ElasticMatrix(Edges,le,ln))));

pause(2);
close all;

display('Constructing principal tree without shrinking grammar');

data = load('test_data/tree23/tree23_inflated.data');
tic;
[NodePositions,Edges,ReportTable] = computeElasticPrincipalGraph(data,50,'GrowGrammars',[{'bisectedge','addnode2node'}],'ShrinkGrammars',[]);
toc;
ln = ones(1,max(max(Edges)));
le = ones(1,size(Edges,1));
display(sprintf('Barcode for the tree = %s',getPrimitiveGraphStructureBarCode(Encode2ElasticMatrix(Edges,le,ln))));

pause(2);
close all; drawnow;

display('Constructing principal tree with mild braching control');

data = load('test_code/branchingTest/thick_turn.data');
nnodes = 50; computeElasticPrincipalGraph(data,nnodes,'Plots',2); 
x = get(gcf,'Position'); 
set(gcf,'Position',[x(1)-x(3) x(2) x(3) x(4)]); 
computeElasticPrincipalGraph(data,nnodes,'BranchingControls',[0.01,1],'Plots',2);

pause(3);
close all; drawnow;

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
lcm = drawPieChartsMetroMap(NodePositions, irisz, pointlabels, NodesMM);


pause(2);
close all;

display('Showing Metro Map Layout, specifying colors for pie-charts');


lcm('Iris-setosa') = [1 1 0]; % Setosa will be yellow
lcm('Iris-versicolor') = [0 1 1]; % Versicolor will be light blue
lcm('Iris-virginica') = [1 0 1]; % Virginica will be purple
close all; 
drawGraph2D(NodesMM,Edges,'ShowClusterNumbers',0); 
drawPieChartsMetroMap(NodePositions, irisz, pointlabels, NodesMM, 'LabelColorMap', lcm);

pause(2);
close all; drawnow;

%%%%%%%%%%%%%%%%% Test 4
display('==================== Test 4 ========================');

display('Preparing data from Human Genome Diversity Project');

hgdptable = importdata('test_data/HGDP_SNP/hgdp_PC3_annotated.txt');
hgdp = hgdptable.data; 
labelRegions = (hgdptable.textdata(2:end,3))';
labelNations = (hgdptable.textdata(2:end,2))';

lcm = containers.Map;
for i=1:size(labelNations,2)
     region = char(labelRegions(i)); label = char(labelNations(i));
 if strcmp(region,'AFRICA') lcm(label) = [0.5 0 0]; end
 if strcmp(region,'SOUTH_CENTRAL_ASIA') lcm(label) = [1 0.5 0]; end
 if strcmp(region,'SOUTH_AMERICA') lcm(label) = [1 0 0]; end
 if strcmp(region,'OCEANIA') lcm(label) = [0 1 1]; end
 if strcmp(region,'EAST_ASIA') lcm(label) = [1 1 0]; end
 if strcmp(region,'EUROPE') lcm(label) = [0 1 0]; end
 if strcmp(region,'NEAR_EAST') lcm(label) = [0.75 0.75 0.75]; end
end
lcm('Sardinian') = [0 1 0]; lcm('Tuscan') = [0.1 1 0.1]; lcm('Italian') = [0.3 1 0.1]; 
lcm('Basque') = [0 0.9 0]; lcm('French') = [0 0.8 0]; lcm('Orcadian') = [0 0.6 0];
lcm('Russian') = [0.1 1 0]; lcm('Adygei') = [0 0.8 0.3];


display('Computing the graph in three stages');
pause(2);

display('First stage: constructing principal curve');
[NodePositions,Edges] = computeElasticPrincipalCurve(hgdp,20,'Plots',2); 
display('Second stage: constructing principal principal graph');
[NodePositions,Edges] = computeElasticPrincipalGraph(hgdp,40,'Plots',2,'Lambda',0.01,'Mu',0.1,'InitGraph',struct('InitNodes',NodePositions,'InitEdges',Edges));
display('Third stage: constructing robust principal principal graph approximating fine details');
[NodePositions,Edges] =  computeElasticPrincipalGraph(hgdp,60,'InitGraph',struct('InitNodes',NodePositions,'InitEdges',Edges),'Lambda',0.0002,'Mu',0.0001,'TrimmingRadius',0.25,'Plots',2);
NodesMM = computeMetroMapLayout(NodePositions,Edges);
figure;
drawGraph2D(NodesMM,Edges,'ShowClusterNumbers',0);
drawPieChartsMetroMap(NodePositions, hgdp, labelNations, NodesMM, 'LabelColorMap', lcm);

pause(2);
close all; drawnow;

display('Showing the piecharts together with tree and data points in projection onto the PCA plane');

[LabelColorMap,partition] = drawPieChartsMetroMap(NodePositions, hgdp, labelNations, NodesMM); 
close all; 
drawPieChartsProc(NodePositions,partition,labelNations,NodesMM,'LabelColorMap',lcm,'ScaleCharts',0.8); 
drawGraph2D(NodePositions,Edges,'NodeSizes',zeros(size(NodePositions,1),1),'ShowClusterNumbers',0); 
for i=1:size(hgdp,1)
    if ~strcmp(labelNations(i),'_')
        plot(hgdp(i,1),hgdp(i,2),'ks','MarkerFaceColor',lcm(char(labelNations(i))),'MarkerSize',6);
        hold on;
    end;
end;

pause(5);
close all; drawnow;
    