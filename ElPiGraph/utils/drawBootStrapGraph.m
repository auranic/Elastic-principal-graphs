function [] = drawBootStrapGraph(AllNodePositions,AllAdjacencyMatrices)

N = size(AllNodePositions,3);

for i=1:N
    a = AllAdjacencyMatrices(:,:,i);
    np = AllNodePositions(:,:,i);
    [row, col] = find(triu(a, 1));
    ed = [row, col];
    drawGraph2D(np,ed,'ShowClusterNumbers',0); hold on;
end

end

