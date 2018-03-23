    close all;
    NodesMM = computeMetroMapLayout(NodePositions,Edges);
    drawGraph2D(NodesMM,Edges,'ShowClusterNumbers',0);
    drawPieChartsMetroMap(NodePositions,hgdp,labelNations,'LabelColorMap',lcm);
