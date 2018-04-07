function drawGraph2D(NodePositions,Edges,varargin)
%drawGraph2D draw graph with nodes coordinates NodePositions and edges
%specified in Edges. 
%
%   IMPORTANT This procedure does not create figure!
%
%Input features
%   NodePositions is k-by-2 array of x and y coordinates of nodes.
%   Edges is k-by-2 matrix of integers. Edges(i,1) and Edges(i,2) specify
%       numbers of two vertex of i-th edge.
%   'Name',Value pairs can customize view:
%       'ShowClusterNumbers' is trigger to write number of each node.
%           Default value is 1 (write)
%       'NodeSizes' preceed array of node sizes. Must have k elements.
%           Default values is array of zeros.
%       'LineWidth' is width of line in graph. Default value is 1;

    % Calculate auxiliary variables
    NodeNum = size(NodePositions,1);
    EdgeNum = size(Edges,1);
    % Specify default values
    NodeSizes = zeros(NodeNum,1);
    LineWidth = 1;
    LineColor = 'r';
    showClusterNumbers = 1;
    % Parse input features
    for i = 1:2:length(varargin)
        if strcmpi(varargin{i}, 'ShowClusterNumbers')
            showClusterNumbers = varargin{i + 1};
        elseif strcmpi(varargin{i}, 'NodeSizes')
            NodeSizes = varargin{i + 1};
        elseif strcmpi(varargin{i}, 'LineWidth')
            LineWidth = varargin{i + 1};
        elseif strcmpi(varargin{i}, 'LineColor')
            LineColor = varargin{i + 1};
        end
    end

    % Check if Edges massif contains zero (means nodes counting starts from
    % zero)
    if any(Edges(:) == 0);
        Edges = Edges + 1;
    end
    
    %Check size of NodeSizes
    if length(NodeSizes) ~= NodeNum
        error(['Parameter NodeSizes must have the same number of',...
            ' elements as number of rows in matrix NodePositions']);
    end
    
    % Rescale node sizes
    ns = double(NodeSizes);
    scale = 20 * sqrt(sum(std(NodePositions(:,1:2)))) / max(ns);
    ns = ns * scale;
    np = double(NodePositions);
    
    % Start drawing. Draw nodes
    hold on;
    for i=1:NodeNum
        if ns(i)>0
            plot(np(i,1),np(i,2),'bo','MarkerSize',...
                ns(i),'LineWidth',LineWidth); hold on;
        else
            plot(np(i,1),np(i,2),'bs','MarkerSize',1,...
                'LineWidth',LineWidth); hold on;
        end
        if showClusterNumbers
            text(np(i,1),np(i,2),sprintf('%i',i),'FontSize',10);
        end
    end
    
    % Draw edges
    for i=1:EdgeNum
        plot([np(Edges(i,1),1) np(Edges(i,2),1)],...
             [np(Edges(i,1),2) np(Edges(i,2),2)],...
             'r-','LineWidth',LineWidth,'Color',LineColor);
    end
    rotate3d off;
end