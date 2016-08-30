LabelColorMap = createLabelColorMap('..\test_data\hgdp_sample_ordered_classes.txt',2);

labels = getColumnOfTable('..\test_data\hgdp_sample_ordered_classes.txt',2);
labelsRegion = getColumnOfTable('..\test_data\hgdp_sample_ordered_classes.txt',3);
data = load('..\test_data\hgdp_PC3.data');

for i=1:size(data,1)
    region = char(labelsRegion(i));
    label = char(labels(i));
    if strcmp(region,'AFRICA')
        LabelColorMap(label) = [0.5 0 0];
    end
    if strcmp(region,'SOUTH_CENTRAL_ASIA')
        LabelColorMap(label) = [1 0.5 0];
    end
    if strcmp(region,'SOUTH_AMERICA')
        LabelColorMap(label) = [1 0 0];
    end
    if strcmp(region,'OCEANIA')
        LabelColorMap(label) = [0 1 1];
    end
    if strcmp(region,'EAST_ASIA')
        LabelColorMap(label) = [1 1 0];
    end
    if strcmp(region,'EUROPE')
        LabelColorMap(label) = [0 1 0];
    end
    if strcmp(region,'NEAR_EAST')
        LabelColorMap(label) = [0.75 0.75 0.75];
    end
        
end

LabelColorMap('Druze') = [0.5 0.5 0.5];
LabelColorMap('Palestinian') = [0.7 0.7 0.7];
LabelColorMap('Bedouin') = [0.8 0.8 0.8];

LabelColorMap('Sardinian') = [0 1 0];
LabelColorMap('Tuscan') = [0.1 1 0.1];
LabelColorMap('Italian') = [0.3 1 0.1];
LabelColorMap('Basque') = [0 0.9 0];
LabelColorMap('French') = [0 0.8 0];
LabelColorMap('Orcadian') = [0 0.6 0];

LabelColorMap('Russian') = [0.1 1 0];
LabelColorMap('Adygei') = [0 0.8 0.3];

LabelColorMap('Kalash') = [1 0.3 0];
LabelColorMap('Pathan') = [1 0.4 0];

LabelColorMap('Balochi') = [1 0.7 0];
LabelColorMap('Brahui') = [0.9 0.6 0];
LabelColorMap('Makrani') = [0.8 0.3 0];
LabelColorMap('Sindhi') = [0.7 0.4 0];
LabelColorMap('Burusho') = [0.8 0.5 0];


LabelColorMap('Mozabite') = [0.5 0 0];

close all;

figure;
drawGraph2D(NodePositions,Edges,'ShowClusterNumbers',0);

for i=1:size(data,1)
    label = char(labels(i));
    color = 'b';
    if ~strcmp(label,'_')
    if LabelColorMap.isKey(label)
        color = LabelColorMap(label);
    end
    plot(data(i,1),data(i,2),'ks','MarkerFaceColor',color,'MarkerSize',6);
    end
end

figure; drawPieCharts(nodesMM,TaxonMap,labels,LabelColorMap,'ScaleCharts',1); drawGraph2D(nodesMM,Edges,'ShowClusterNumbers',0); 
%xlim([-2 -0.5]); ylim([-0.35 0.7]); 
axis off;