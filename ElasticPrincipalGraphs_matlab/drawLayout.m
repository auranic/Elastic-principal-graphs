function drawLayout(nodes,nodes2d,edges)

figure; drawGraph2D(nodes2d,edges); figure; [v,u] = pca(nodes); drawGraph2D(u(:,1:2),edges);