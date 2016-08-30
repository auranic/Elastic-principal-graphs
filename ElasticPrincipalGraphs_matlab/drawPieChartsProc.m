function [LabelColorMap] = drawPieChartsProc(NodePositions, TaxonMap, Labels, varargin)

makeColorMap = 1;

    for i=1:length(varargin)
        if strcmpi(varargin{i},'LabelColorMap')
            makeColorMap = 0;
            lcm = varargin(i+1);
            LabelColorMap = lcm{1};
        end
    end

if makeColorMap
    LabelColorMap = createLabelColorMapList(Labels);    
else
    tempMap = createLabelColorMapList(Labels);
    tempLabels = LabelColorMap.keys;
    for i=1:size(tempLabels,2)
        s = char(tempLabels(i));
        tempMap(s) = LabelColorMap(s);
    end
    LabelColorMap = tempMap;
end

scale = sqrt(sum(std(NodePositions(:,1:2))));
%NodeNum = size(NodePositions,1);
NodeNum = TaxonMap.Count;
NodeSizes = zeros(NodeNum,1);
NodeLabels = TaxonMap.keys;
for i=1:size(NodeLabels,2)
    k = str2num(char(NodeLabels(i)));
    NodeSizes(k) = size(TaxonMap(char(NodeLabels(i))),1);
end

%ns = cast(NodeSizes,'DOUBLE');
%ns = NodeSizes;
[ns,ind] = sort(NodeSizes,'descend');
%ns = NodeSizes;
%ind = 1:NodeNum;

scaleCharts = 1;

    for i=1:length(varargin)
        if strcmpi(varargin{i},'ScaleCharts')
            scaleCharts = varargin{i+1};
        elseif strcmpi(varargin{i},'')
            something = varargin{i+1};
        end
    end


maxNS = max(ns);
ns = ns/maxNS*scale/5*scaleCharts;


for i=1:size(NodePositions,1)
    
    taxon = TaxonMap(num2str(ind(i)));
    countsMap = containers.Map;
    
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
    
    propLabels = keys(countsMap);
    numLabels = size(propLabels,2);
    props = zeros(1,numLabels);
    colors = zeros(numLabels,3);
    
    for j=1:numLabels
        label = char(propLabels(j));
        props(j) = countsMap(label);
        colors(j,:) = LabelColorMap(label);
    end
    
    drawPieChart(NodePositions(ind(i),1),NodePositions(ind(i),2),ns(i),props,colors); hold on;
    
end

axis equal; axis auto;