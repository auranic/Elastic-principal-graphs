% In this example we assume the number of trajectories known (3)
% It is relatively simple to make an algorithm which will find the number
% of trajectories automatically, by checking if the set of captured points
% is empty.
% The trimming radius can be also determined automatically.
% In case of 'thick' trajectories and non-orthogocal trajectory
% intersections, the performance of the algorithm decreases.
% There are various ideas for improving it: for example, 
% in the areas of the captured points, the data can be downsampled,
% in order to decrease density (but not remove the data points completely,
% one will need 'bridges' between trajectory fragments.

TrimmingRadius = 4;
Lambda = 0.00001;
Mu = 0.2;
Nnodes = 200;

%%%%%%%%%% Reading the image data and converting to a cloud of points in 2D
[A,map] = imread('tm3_tiny.png'); 
A1 = 255-A; 
[x,y] = find(A1(:,:,1)>0); 
plot(y,x,'k.'); axis equal; tm = [y x];
%break;

%%%%%%%%%%%%%%%%%%%%%%
% First trajectory

tic;
[np,ed] = computeRobustElasticPrincipalGraph(tm,Nnodes,TrimmingRadius,'Plots',0,'Lambda',Lambda,'Mu',Mu,...
    'GrowGrammars',[{'bisectedge';'addnode2terminalnode'}],'ShrinkGrammars',[]); 
toc; 

figure; 
drawSimple2DIllustration(tm,np,ed,'Color','r'); 

SquaredX = sum(tm.^2, 2);
MaxBlockSize = 10000;

%%%%%%%%%%%%%%%%%%%%%%
% Second trajectory

% initialization from the points which are not yet captured
[partition] = PartitionData(tm, np, MaxBlockSize,SquaredX, TrimmingRadius*3);
inds = partition==0;
[np0,ed0] = computeRobustElasticPrincipalGraph(tm(inds,:),2,TrimmingRadius,'Plots',0,'Lambda',Lambda,'Mu',Mu); 
figure;
plot(tm(inds,1),tm(inds,2),'kx'); drawGraph2D(np0,ed0,'LineColor','m','LineWidth',10,'ShowClusterNumbers',0);  drawnow;

% compute second trajectory
tic; [np1,ed1] = computeElasticPrincipalCurve(tm,Nnodes,'Plots',0,'TrimmingRadius',TrimmingRadius,'Lambda',Lambda,'Mu',Mu,...
    'InitGraph',struct('InitNodes',np0,'InitEdges',ed0),...
    'GrowGrammars',[{'bisectedge';'addnode2terminalnode'}],'ShrinkGrammars',[]); toc;
    
hold on; drawGraph2D(np1,ed1,'LineColor','b','LineWidth',3,'ShowClusterNumbers',0); 

%%%%%%%%%%%%%%%%%%%%%%
% Third trajectory

% initialization from the points which are not yet captured
temp = [np;np1];
[partition] = PartitionData(tm, temp, MaxBlockSize,SquaredX, TrimmingRadius*3);
inds = partition==0;
[np0,ed0] = computeRobustElasticPrincipalGraph(tm(inds,:),2,TrimmingRadius,'Plots',0,'Lambda',Lambda,'Mu',Mu); 
figure;
plot(tm(inds,1),tm(inds,2),'kx'); drawGraph2D(np0,ed0,'LineColor','m','LineWidth',10,'ShowClusterNumbers',0);  drawnow;

% compute third trajectory
tic; [np2,ed2] = computeElasticPrincipalCurve(tm,Nnodes,'Plots',0,'TrimmingRadius',TrimmingRadius,'Lambda',Lambda,'Mu',Mu,...
    'InitGraph',struct('InitNodes',np0,'InitEdges',ed0),...
    'GrowGrammars',[{'bisectedge';'addnode2terminalnode'}],'ShrinkGrammars',[]); toc;


% draw final image
figure;
drawSimple2DIllustration(tm,np,ed,'Color','r'); 
hold on; drawGraph2D(np1,ed1,'LineColor','b','LineWidth',3,'ShowClusterNumbers',0); 
hold on; drawGraph2D(np2,ed2,'LineColor','g','LineWidth',3,'ShowClusterNumbers',0); 

% cluster data by trajectories
[partition1] = PartitionData(tm, np, MaxBlockSize,SquaredX, TrimmingRadius);
[partition2] = PartitionData(tm, np1, MaxBlockSize,SquaredX, TrimmingRadius);
[partition3] = PartitionData(tm, np2, MaxBlockSize,SquaredX, TrimmingRadius);
ind1 = partition1>0;
ind2 = partition2>0;
ind3 = partition3>0;
figure;
plot(tm(ind1,1),tm(ind1,2),'ro'); hold on;
plot(tm(ind2,1),tm(ind2,2),'bo'); hold on;
plot(tm(ind3,1),tm(ind3,2),'go'); hold on;
axis off; axis equal;