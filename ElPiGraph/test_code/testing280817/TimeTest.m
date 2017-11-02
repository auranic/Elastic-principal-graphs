nRep = 1;
tic;
for k = 1:nRep
    [EmbeddedNodePositions, ElasticEnergy, partition, MSE,EP,RP] =...
        PrimitiveElasticGraphEmbedment(X, NodePositions, ElasticMatrix, 'TrimmingRadius', 15, 'verbose', 1);
end
disp(toc/nRep);
