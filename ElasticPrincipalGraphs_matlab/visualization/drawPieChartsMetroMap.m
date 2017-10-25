function [LabelColorMap, graph, TaxonMap] =...
    drawPieChartsMetroMap(NodePositions, Edges, data, labels, varargin)
%drawPieChartsMetroMap draw metro map with pie cart in each node
%
%   IMPORTANT This procedure does not create figure!
%
%Input features
%   NodePositions is k-by-2 array of x and y coordinates of nodes.
%   Edges is k-by-2 matrix of integers. Edges(i,1) and Edges(i,2) specify
%       numbers of two vertex of i-th edge.
%   data is n-by-m matrix data (each row of matrix contains one observation).
%   labels is vector of labels (one label for each observation).
%   'Name',Value pairs can customize view:
%       'LabelColorMap' is list of labels for color map. Default value is
%           automative creation of colour map in accordance with Labels
%       'ScaleCharts' is number to scale charts
%   

    % label must be row vector of labels for all observations
    if size(labels,1) > 1
        labels = labels';
    end
    if size(labels,1) > 1 || size(labels,2) ~= size(data,1)
        error(['labels must be a vector with the same number of',...
            ' elements as number of observations (rows) in data']);
    end
    % Calculate positions of nodes in metro map
    nodesMM = computeMetroMapLayout(NodePositions,Edges);
    % Calculate graph
    graph = makeGraph(NodePositions,Edges);
    % Form taxon map
    TaxonMap = getTaxonMap(graph,data);
    % Draw diagram
    hold on;
    [LabelColorMap] =...
        drawPieChartsProc(nodesMM, TaxonMap, labels, varargin{:});
end