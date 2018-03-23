function drawMetroMap(NodePositions, Edges, varargin)
%drawMetroMap draw metro map by NodePositions and Edges
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
%       'DrawPieChart' is true to draw pie chart and false otherwise.
%       'Labels' is is vector of labels (one label for each observation).
%           Must be preseted to draw pie chart 
%       'Data' is is n-by-m matrix data (each row of matrix contains one
%           observation). Must be preseted to draw pie chart 
%       'LabelColorMap' is list of labels for color map. Default value is
%           automative creation of colour map in accordance with Labels.
%           Must be preseted to draw pie chart.
%       'ScaleCharts' is number to scale charts
%


    % Specify default values
    NodeSizes = zeros(size(NodePositions,2),1);
    LineWidth = 1;
    showClusterNumbers = 1;
    drawPieChart = 0;
    data = 0;
    labels = 0;
    % Parse input features
    for i = 1:2:length(varargin)
        if strcmpi(varargin{i}, 'ShowClusterNumbers')
            showClusterNumbers = varargin{i + 1};
        elseif strcmpi(varargin{i}, 'NodeSizes')
            NodeSizes = varargin{i + 1};
        elseif strcmpi(varargin{i}, 'LineWidth')
            LineWidth = varargin{i + 1};
        elseif strcmpi(varargin{i}, 'DrawPieChart')
            drawPieChart = varargin{i + 1};
        elseif strcmpi(varargin{i}, 'Data')
            data = varargin{i + 1};
        elseif strcmpi(varargin{i}, 'Labels')
            labels = varargin{i + 1};
        end
    end
    
    % Calculate layout of metro map
    NodesMM = computeMetroMapLayout(NodePositions,Edges);
    %Create figure
    figure;
    hold on;
    drawGraph2D(NodesMM,Edges,'ShowClusterNumbers', showClusterNumbers,...
        'NodeSizes', NodeSizes ,'LineWidth', LineWidth);
    if drawPieChart
        if isscalar(data) || isscalar(labels)
            error('Data and labels must be presented to draw pie chart');
        end
        drawPieChartsMetroMap(NodePositions, data, labels, varargin{:});
   end
   title('Metro map layout of the principal tree','FontSize',20);
   hold off;
end