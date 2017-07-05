function [LabelColorMap] =...
    drawPieChartsProc(NodePositions, TaxonMap, Labels, varargin)
%drawPieChartsProc is procedure to draw pie chats in each nodeof graph
%
%   IMPORTANT This procedure does not create figure!
%
%Input features:
%   NodePositions is k-by-2 array of x and y coordinates of nodes.
%   TaxonMap is map of data points. Each element of map contains list of
%       data points associated with this node 
%   Labels is vector of labels (one label for each observation).
%   'Name',Value pairs can customize view:
%       'LabelColorMap' is list of labels for color map. Default value is
%           automative creation of colour map in accordance with Labels
%       'ScaleCharts' is number to scale charts

    % Parse optional parameters
    makeColorMap = 1;
    scaleCharts = 1;
    for i=1:2:length(varargin)
        if strcmpi(varargin{i},'LabelColorMap')
            makeColorMap = 0;
            lcm = varargin(i+1);
            LabelColorMap = lcm{1};
        elseif strcmpi(varargin{i},'ScaleCharts')
            scaleCharts = varargin{i+1};
        end
    end
    
    % Creater label color map in accordance with user specification
    if makeColorMap
        % Automative creation
        LabelColorMap = createLabelColorMapList(Labels);    
    else
        % Use map created by user ???
        tempMap = createLabelColorMapList(Labels);
        tempLabels = LabelColorMap.keys;
        for i=1:size(tempLabels,2)
            s = char(tempLabels(i));
            tempMap(s) = LabelColorMap(s);
        end
        LabelColorMap = tempMap;
    end
    % Calculate scale for charts
    scale = sqrt(sum(std(NodePositions(:,1:2))));
    % Calculate node sizes from TaxonMap
    NodeNum = TaxonMap.Count;
    NodeSizes = zeros(NodeNum,1);
    NodeLabels = TaxonMap.keys;
    for i=1:size(NodeLabels,2)
        k = str2num(char(NodeLabels(i)));
        NodeSizes(k) = size(TaxonMap(char(NodeLabels(i))),1);
    end
    % Sort sizes of nodes
    [ns,ind] = sort(NodeSizes,'descend');
    % Rescale data
    ns = ns / max(ns) * scale / 5 * scaleCharts;
    %Draw pie charts
    for i=1:size(NodePositions,1)
        %Get taxon
        taxon = TaxonMap(num2str(ind(i)));
        countsMap = containers.Map;
        % Calculate number of cases of each class
        for j=1:size(taxon,1)
            label = char(Labels(taxon(j)));
            if ~strcmp(label,'_')
            if ~countsMap.isKey(label)
                countsMap(label) = 0;
            end
            count = countsMap(label)+1;
            countsMap(label) = count;
            end
        end
        % Get proportion label and colour for each class
        propLabels = keys(countsMap);
        numLabels = size(propLabels,2);
        props = zeros(1,numLabels);
        colors = zeros(numLabels,3);
        for j=1:numLabels
            label = char(propLabels(j));
            props(j) = countsMap(label);
            colors(j,:) = LabelColorMap(label);
        end
        %Draw chart
        drawPieChart(NodePositions(ind(i),1), NodePositions(ind(i),2),...
            ns(i),props,colors); 
    end
    % Andrei, I am not sure that it is good idea to put this actions in
    % this function ???
    axis equal; axis auto;
end