function PCAView( Nodes, Edges, data, pc1, pc2, pc1FVE, pc2FVE )
%PCAView draw the dots from data and graph in the space of the two selected
%principal components. Data will be centralized before projection. Mean
%vector of data will be subtracted from nodes positions.
%
%Input features
%   Nodes is array of nodes positions
%   Edges is k-by-2 matrix of integers. Edges(i,1) and Edges(i,2) specify
%       numbers of two vertex of i-th edge.
%   data is n-by-m matrix data (each row of matrix contains one observation).
%   pc1 is optional argument: 
%       If it is omitted then 1 is used. 
%       If it is integer number then it is number of principal component to
%           use.
%       If it is m-by-1 array then it is vector to calculate projections
%   pc2 is optional argument: 
%       If it is omitted and pc1 is integer number then pc1 + 1 is used. 
%       If it is omitted and pc1 is vector then 1 is used. 
%       If it is integer number then it is number of principal component to
%           use.
%       If it is m-by-1 array then it is vector to calculate projections
%   pc1FVE is optional argument. Used as part of axis title if specified.
%       If pc1 is number of component then this value is calculated.
%   pc2FVE is optional argument. Used as part of axis title if specified.
%       If pc2 is number of component then this value is calculated.

    % Check input arguments
    if nargin < 3
        error('At least Nodes, Edges, data must be specified');
    end
    % Centralise data
    means = mean(data);
    data = bsxfun(@minus, data, means);
    Nodes = bsxfun(@minus, Nodes, means);
    % Check pc1 and pc2
    if nargin < 5
        pc1 = 1;
        pc2 = 2;
        pc1FVE = 0;
        pc2FVE = 0;
    elseif nargin < 6
        if length(pc1) == 1
            pc2 = pc1 + 1;
            pc2FVE = 0;
        else
            pc2 = 1;
        end
    elseif nargin < 7
        pc1FVE = 0;
        pc2FVE = 0;
    elseif nargin < 8
        pc2FVE = 0;
    end
    
    % Calculate principal components if necessary
    if length(pc1) == 1 || length(pc2) == 1
        [vglobal, uglobal, explainedVariances] = pca(data);
    end
    % Form coordinates to draw
    if length(pc1) == 1
        xData = uglobal(:,pc1);
        xNodes = Nodes*vglobal(:,pc1);
        pc1FVE = explainedVariances(pc1)/sum(explainedVariances);
    else
        xData = data*pc1(:);
        xNodes = Nodes*pc1(:);
    end
    if length(pc2) == 1
        yData = uglobal(:,pc2);
        yNodes = Nodes*vglobal(:,pc2);
        pc1FVE = explainedVariances(pc2)/sum(explainedVariances);
    else
        yData = data*pc2(:);
        yNodes = Nodes*pc2(:);
    end
    % Get size of nodes
    SquaredX = sum(data.^2, 2);
    MaxBlockSize = 10000;
    [partition] = PartitionData(data,Nodes,MaxBlockSize,SquaredX);
    NodeSizes = ones(size(Nodes,1),1);
     for i=1:size(Nodes,1)
         NodeSizes(i) = sum(partition==i)+1;
     end
    
    % Create figure
    figure;
    
    % Draw data
    %plot(xData,yData,'ko','MarkerSize',2);
    % Draw data clustered by branches
    em = MakeUniformElasticMatrix(Edges,1,1);
    [node_partition, internal_flag, star_centers] = partition_nodes_by_branch(em);
    inds1 = find(node_partition==0);
    inds = [];
    for i=1:size(inds1,1)
        inds = [inds;find(partition==inds1(i))];
    end
    plot(xData(inds),yData(inds),'ko','MarkerSize',2); hold on;
    for i=1:max(node_partition)
        labels(i) = {int2str(i)};
    end
    [LabelColorMap] = createLabelColorMapList(labels);
    
    for i=1:max(node_partition)
        inds1 = find(node_partition==i);
        inds = [];
        for j=1:size(inds1,1)
            inds = [inds;find(partition==inds1(j))];
        end
        color = LabelColorMap(char(int2str(i)));
        plot(xData(inds),yData(inds),'ko','MarkerSize',2,'MarkerFaceColor',color,'MarkerEdgeColor',color);
    end
    
    
    % Draw tree with specified sizes of nodes
    drawGraph2D([xNodes,yNodes],Edges,'NodeSizes',NodeSizes);
    
    
    title('PCA view of principal graph','FontSize',20);
    if pc1FVE>0
        xlabel(sprintf('PCx(%2.2f%%)',pc1FVE*100),'FontSize',20);
    else
        xlabel('PCx','FontSize',20);
    end
    if pc2FVE>0
        ylabel(sprintf('PCy(%2.2f%%)',pc2FVE*100),'FontSize',20);
    else
        ylabel('PCy','FontSize',20);
    end
end

