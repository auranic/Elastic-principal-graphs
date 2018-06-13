iristable = importdata('test_data\iris\iris.txt'); 
iris = iristable.data;
irisz = zscore(iris);
pointlabels = (iristable.textdata(2:end,1))';

[v,u,s] = pca(irisz);

lcm = containers.Map;
lcm('Iris-setosa') = [0 0 1]; % Setosa will be yellow
lcm('Iris-versicolor') = [0 1 0]; % Versicolor will be light blue
lcm('Iris-virginica') = [1 0 1]; % Virginica will be purple

mu = 0.2;

subplot(2,2,1);

[np,ed] = computeElasticPrincipalGraph(irisz,30,'BranchingControls',[0 1],'Plots',0,'Mu',mu);
npu = np*v;


plot(u(1:49,1),u(1:49,2),'bo'); hold on;
plot(u(50:99,1),u(50:99,2),'ko','MarkerEdgeColor',[0 0.7 0]); hold on;
plot(u(100:149,1),u(100:149,2),'ko','MarkerEdgeColor',[1 0 1]); hold on;

xlabel('PC1'); ylabel('PC2'); title('Alpha = 0');
drawGraph2D(npu(:,1:2),ed,'ShowClusterNumbers',0,'LineWidth',5); hold on;
[LabelColorMap,partition] = drawPieChartsMetroMap(np, irisz, pointlabels, npu,'LabelColorMap',lcm); 


subplot(2,2,2);

[np,ed] = computeElasticPrincipalGraph(irisz,30,'BranchingControls',[0.005 1],'Plots',0,'Mu',mu);
npu = np*v;


plot(u(1:49,1),u(1:49,2),'bo'); hold on;
plot(u(50:99,1),u(50:99,2),'ko','MarkerEdgeColor',[0 0.7 0]); hold on;
plot(u(100:149,1),u(100:149,2),'ko','MarkerEdgeColor',[1 0 1]); hold on;

xlabel('PC1'); ylabel('PC2'); title('Alpha = 0.005');
drawGraph2D(npu(:,1:2),ed,'ShowClusterNumbers',0,'LineWidth',5); hold on;
[LabelColorMap,partition] = drawPieChartsMetroMap(np, irisz, pointlabels, npu,'LabelColorMap',lcm); 


subplot(2,2,3);

[np,ed] = computeElasticPrincipalGraph(irisz,30,'BranchingControls',[0.01 1],'Plots',0,'Mu',mu);
npu = np*v;


plot(u(1:49,1),u(1:49,2),'bo'); hold on;
plot(u(50:99,1),u(50:99,2),'ko','MarkerEdgeColor',[0 0.7 0]); hold on;
plot(u(100:149,1),u(100:149,2),'ko','MarkerEdgeColor',[1 0 1]); hold on;

xlabel('PC1'); ylabel('PC2'); title('Alpha = 0.01');
drawGraph2D(npu(:,1:2),ed,'ShowClusterNumbers',0,'LineWidth',5); hold on;
[LabelColorMap,partition] = drawPieChartsMetroMap(np, irisz, pointlabels, npu,'LabelColorMap',lcm); 

subplot(2,2,4);

[np,ed] = computeElasticPrincipalGraph(irisz,30,'BranchingControls',[1 1],'Plots',0,'Mu',mu);
npu = np*v;


plot(u(1:49,1),u(1:49,2),'bo'); hold on;
plot(u(50:99,1),u(50:99,2),'ko','MarkerEdgeColor',[0 0.7 0]); hold on;
plot(u(100:149,1),u(100:149,2),'ko','MarkerEdgeColor',[1 0 1]); hold on;

xlabel('PC1'); ylabel('PC2'); title('Alpha = 1');
drawGraph2D(npu(:,1:2),ed,'ShowClusterNumbers',0,'LineWidth',5); hold on;
[LabelColorMap,partition] = drawPieChartsMetroMap(np, irisz, pointlabels, npu,'LabelColorMap',lcm); 

