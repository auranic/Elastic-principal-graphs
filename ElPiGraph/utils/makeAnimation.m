function [NodePositions,Edges] = makeAnimation(data,numNodes,filename,varargin)

Lambda = 0.01;
Mu = 0.1;

    for i=1:2:length(varargin)
        if strcmpi(varargin{i}, 'Lambda')
            Lambda = varargin{i + 1}; 	
        elseif strcmpi(varargin{i}, 'Mu')
            Mu = varargin{i + 1}; 	
        end
    end


[NodePositions,Edges,ReportTable] = ElPrincGraph(data,4,Lambda,Mu,varargin{:});
plot(data(:,1),data(:,2),'ko'); hold on; drawGraph2D(NodePositions(:,1:2),Edges,'NodeSizes',ones(size(NodePositions,1),1),'ShowClusterNumbers',0,'LineWidth',3);  axis off;
drawnow;
saveas(gcf,sprintf('%s%2i.png',filename,2));

k = size(varargin,2);
varargin{k+1} = 'InitNodePositions';
varargin{k+2} = NodePositions;
varargin{k+3} = 'InitElasticMatrix';
varargin{k+4} = MakeUniformElasticMatrix(Edges,Lambda,Mu);

for i=size(NodePositions,1)+1:numNodes
    varargin{k+2} = NodePositions;
    varargin{k+3} = Edges;
    [NodePositions,Edges,ReportTable] = ElPrincGraph(data,i,Lambda,Mu,varargin{:});
    hold off; plot(data(:,1),data(:,2),'ko'); hold on; 
    drawGraph2D(NodePositions(:,1:2),Edges,'NodeSizes',ones(size(NodePositions,1),1),'ShowClusterNumbers',0,'LineWidth',3); axis off;
    drawnow;
    saveas(gcf,sprintf('%s%2i.png',filename,i));
end


end

