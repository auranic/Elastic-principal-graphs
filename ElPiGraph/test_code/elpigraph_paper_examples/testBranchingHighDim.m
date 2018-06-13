clear all;

tab = importdata('dim10.tsv');
data = tab.data(:,1:10);
dataz = zscore(data);

nnodes = 40;

branches = tab.data(:,11);
[v,u,s] = pca(dataz);

inds(1,:) = branches==1;
inds(2,:) = branches==2;
inds(3,:) = branches==3;
inds(4,:) = branches==4;
inds(5,:) = branches==5;
col = zeros(5,3);
col(1,:) = [1;0.5;0.5];
col(2,:) = [0.5;0.5;0.5];
col(3,:) = [0.1;1;0.1];
col(4,:) = [0.5;0.5;1];
col(5,:) = [1;0.5;1];
for i=1:5
plot(u(inds(i,:),1),u(inds(i,:),2),'ko','Color',col(i,:)); hold on;
end

[NodePositions2D,Edges2D] = computeElasticPrincipalGraph(u(:,1:2),nnodes,'BranchingControls',[0.01 1],'Plots',0);
drawGraph2D(NodePositions2D,Edges2D,'LineWidth',6,'ShowClusterNumbers',0);

axis off; axis equal;

figure;

[NodePositions10D,Edges] = computeElasticPrincipalGraph(dataz,nnodes,'BranchingControls',[0.01 1],'Plots',0);
npu = NodePositions10D*v;
for i=1:5
plot(u(inds(i,:),1),u(inds(i,:),2),'ko','Color',col(i,:)); hold on;
end
drawGraph2D(npu(:,1:2),Edges,'LineWidth',6,'ShowClusterNumbers',0);

brmap = containers.Map;

for i=1:5
    brmap(int2str(i)) = col(i,:);
end

clear branches_string;
for i=1:length(branches)
    if(branches(i)==1) branches_string(i)='1'; end;
    if(branches(i)==2) branches_string(i)='2'; end;
    if(branches(i)==3) branches_string(i)='3'; end;
    if(branches(i)==4) branches_string(i)='4'; end;
    if(branches(i)==5) branches_string(i)='5'; end;
end

NodesMM = computeMetroMapLayout(NodePositions10D,Edges);
figure;
drawGraph2D(NodesMM,Edges,'ShowClusterNumbers',0);
drawPieChartsMetroMap(NodePositions10D, dataz, branches_string, NodesMM, 'LabelColorMap', brmap);

axis off; axis equal;

NodesMM = computeMetroMapLayout(NodePositions2D,Edges2D);
figure;
drawGraph2D(NodesMM,Edges2D,'ShowClusterNumbers',0);
drawPieChartsMetroMap(NodePositions2D, u(:,1:2), branches_string, NodesMM, 'LabelColorMap', brmap);

axis off; axis equal;

	