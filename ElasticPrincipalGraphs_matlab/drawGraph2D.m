function drawGraph2D(NodePositions,Edges,varargin)

NodeNum = size(NodePositions,1);
EdgeNum = size(Edges,1);
scale = sqrt(sum(std(NodePositions(:,1:2))));
NodeSizes = zeros(NodeNum,1);


showClusterNumbers = 1;

    for i=1:length(varargin)
        if strcmpi(varargin{i},'ShowClusterNumbers')
            showClusterNumbers = varargin{i+1};
        elseif strcmpi(varargin{i},'NodeSizes')
            NodeSizes = varargin{i+1};
        end
    end


% Check if Edges massif contains zero (means counting nodes starts from
% zero)
tmp = find(Edges==0);
if(size(tmp,1)>0)
    Edges=Edges+1;
end

ns = cast(NodeSizes,'DOUBLE');

maxNS = max(ns);

ns = ns/maxNS*scale*20;

%*scale/10;

np = cast(NodePositions,'DOUBLE');

for i=1:NodeNum
    if ns(i)>0
        plot(np(i,1),np(i,2),'bo','MarkerSize',ns(i)); hold on;
    else
        plot(np(i,1),np(i,2),'bs','MarkerSize',1); hold on;
    end
    if showClusterNumbers
    text(np(i,1),np(i,2),sprintf('%i',i));
    end
end

for i=1:EdgeNum
    plot([NodePositions(Edges(i,1),1) NodePositions(Edges(i,2),1)],[NodePositions(Edges(i,1),2) NodePositions(Edges(i,2),2)],'r-','LineWidth',1);
end

