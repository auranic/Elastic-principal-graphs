edges = [1,2;1,3;1,4;4,5;5,6;5,7];
lambdas = repmat(0.01,size(edges,1),1);
mus = [10,0,0,10,10,0,0];

[ElasticMatrix] = Encode2ElasticMatrix(edges,lambdas,mus);

%ElasticMatrix

[e,l,m] = DecodeElasticMatrix(ElasticMatrix);

%e
%l
%m

X = load('./test_data/tree23/tree23.data');
NNodes = max(max(edges));
NodePositions = X(randperm(size(X,1),NNodes),:);
[NodePositions] = PrimitiveElasticGraphEmbedment(X, NodePositions, ElasticMatrix,'MaxNumberOfIterations',100,'verbose',0);
[ed] = DecodeElasticMatrix(ElasticMatrix);
figure; plot(X(:,1),X(:,2),'b.'); hold on; drawGraph2D(NodePositions,ed);
%[NodePositionArray, ElasticMatrices] = GraphGrammarOperation(NodePositions,ElasticMatrix,X,'removenode');
%[NodePositionArray, ElasticMatrices] = GraphGrammarOperation(NodePositions,ElasticMatrix,X,'shrinkedge');
%[NodePositionArray, ElasticMatrices] = GraphGrammarOperation(NodePositions,ElasticMatrix,X,'bisectedge');
%[NodePositionArray, ElasticMatrices] = GraphGrammarOperation(NodePositions,ElasticMatrix,X,'addnode2node');


% for i=1:size(NodePositionArray,3)
%     [np,energy] = PrimitiveElasticGraphEmbedment(X, NodePositionArray(:,:,i), ElasticMatrices(:,:,i),'MaxNumberOfIterations',100);
%     [ed] = DecodeElasticMatrix(ElasticMatrices(:,:,i));
%     figure; plot(X(:,1),X(:,2),'k.'); hold on; drawGraph2D(np,ed);
%     display(sprintf('Graph %i, Energy=%f',i,energy));
% end



