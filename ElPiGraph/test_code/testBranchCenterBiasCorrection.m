x = load('C:\MyPrograms\_github\Elastic-principal-graphs\ElPiGraph\test_data\LLE_Pinello\AF.X.csv');
n = load('C:\MyPrograms\_github\Elastic-principal-graphs\ElPiGraph\test_data\LLE_Pinello\AF.cluster_centers_.csv');
e = load('C:\MyPrograms\_github\Elastic-principal-graphs\ElPiGraph\test_data\LLE_Pinello\AF.cluster_edges_.csv');

dims = [1:3];
nnodes = size(n,1);
Lambda = 0.01;
Mu = 0.1;
TrimmingRadius = 0.02;

% branching control parameter
Alpha = 0.01;
% factor of softening the stars with k>2 and their edges
beta = 10;

%[Nodes,Edges] = computeElasticPrincipalGraph(x(:,dims),nnodes+30,'BranchingControls',[0.01 1],'Lambda',Lambda,'Mu',Mu,'TrimmingRadius',TrimmingRadius,'Plots',2,'InitGraph',struct('InitNodes',n(:,dims),'InitEdges',e));
%[Nodes,Edges] = computeElasticPrincipalGraph(x(:,dims),nnodes+30,'BranchingControls',[0.0 10],'Lambda',Lambda,'Mu',Mu,'TrimmingRadius',TrimmingRadius,'Plots',2);

Nodes = n(:,dims);
Edges = e;

em = MakeUniformElasticMatrix(Edges, Lambda, Mu);
  L = em - diag(diag(em));
  M = diag(em);
  connectivity = sum(L>0);
  stars = find(connectivity>2);
  M(stars) = Mu/beta;
  for i=1:length(stars)
    leaves = find(L(:,stars(i))>0);
    L(stars(i),leaves) = Lambda/beta;
    L(leaves,stars(i)) = Lambda/beta;
  end
  em1 = L+diag(M);

GraphStruct = struct('InitNodes',Nodes,'InitEdges',Edges,'ElasticMatrix',em1);
[Nodes,Edges] = computeElasticPrincipalGraph(x(:,dims),nnodes+23,'BranchingControls',[Alpha 1],'Lambda',Lambda,'Mu',Mu,'TrimmingRadius',TrimmingRadius,'Plots',2,'InitGraph',GraphStruct);

figure;
% plot(x(:,2),x(:,3),'k.'); hold on; plot(Nodes(:,2),Nodes(:,3),'ro');
% for i=1:size(Edges,1)
%     plot([Nodes(Edges(i,1),2) Nodes(Edges(i,2),2)],[Nodes(Edges(i,1),3) Nodes(Edges(i,2),3)],'b-');
% end

plot3(x(:,1),x(:,2),x(:,3),'k.'); hold on; plot3(Nodes(:,1),Nodes(:,2),Nodes(:,3),'ro');
for i=1:size(Edges,1)
    plot3([Nodes(Edges(i,1),1) Nodes(Edges(i,2),1)],[Nodes(Edges(i,1),2) Nodes(Edges(i,2),2)],[Nodes(Edges(i,1),3) Nodes(Edges(i,2),3)],'b-');
end
