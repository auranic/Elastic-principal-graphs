function makeAnimation(data,numNodes,varargin)


[NodePositions,Edges,ReportTable,cpg] = computeElPT(data,2,varargin{:});
plot(data(:,1),data(:,2),'ko'); hold on; drawGraph2D(NodePositions(:,1:2),Edges,'NodeSizes',ones(size(NodePositions,1)),'ShowClusterNumbers',0,'LineWidth',3);  axis off;
drawnow;
saveas(gcf,sprintf('animation/node%2i.png',2));

k = size(varargin,2);
varargin{k+1} = 'InitGraph';

for i=size(NodePositions,1)+1:numNodes
    varargin{k+2} = NodePositions;
    varargin{k+3} = Edges;
    [NodePositions,Edges,ReportTable,cpg] = computeElPT(data,i,varargin{:});
    hold off; plot(data(:,1),data(:,2),'ko'); hold on; drawGraph2D(NodePositions(:,1:2),Edges,'NodeSizes',ones(size(NodePositions,1)),'ShowClusterNumbers',0,'LineWidth',3); axis off;
    drawnow;
    saveas(gcf,sprintf('animation/node%2i.png',i));
end


end

