cc = load('loop_branch.data');
[np0,ed0] = computeElasticPrincipalCircle(cc,6,'Plots',0);
[np,ed] = computeElasticPrincipalGraph(cc,30,'InitGraph',struct('InitNodes',np0,'InitEdges',ed0),'BranchingControls',[0.01,1],'Plots',2);
