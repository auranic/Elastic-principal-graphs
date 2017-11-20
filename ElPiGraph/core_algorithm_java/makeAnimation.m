function makeAnimation(data, numNodes, filename, varargin)

[NodePositions,Edges,ReportTable] = computeElPT(data,2,varargin{:});
plot(data(:,1),data(:,2),'ko'); hold on; drawGraph2D(NodePositions(:,1:2),Edges,'NodeSizes',ones(size(NodePositions,1),1),'ShowClusterNumbers',0,'LineWidth',3);  axis off;
drawnow;
saveas(gcf,sprintf('%s%2i.png',filename,2));

k = size(varargin,2);
varargin{k+1} = 'InitGraph';

for i=size(NodePositions,1)+1:numNodes
    varargin{k+2} = NodePositions;
    varargin{k+3} = Edges;
    [NodePositions,Edges,ReportTable] = computeElPT(data,i,varargin{:});
    hold off; plot(data(:,1),data(:,2),'ko'); hold on; drawGraph2D(NodePositions(:,1:2),Edges,'NodeSizes',ones(size(NodePositions,1),1),'ShowClusterNumbers',0,'LineWidth',3); axis off;
    drawnow;
    saveas(gcf,sprintf('%s%2i.png',filename,i));
end


end

