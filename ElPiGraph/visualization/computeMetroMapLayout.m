function [NodePositions2D,mml] = computeMetroMapLayout(NodePositions,Edges)

javaclasspath({'.\core_algorithm_java\VDAOEngine.jar'});

mml = vdaoengine.analysis.grammars.MetroMapLayout;

tree = vdaoengine.analysis.grammars.Tree;

mml.principalTree = tree;

e = Edges;

if sum(sum(Edges==0))==0
    e = e-1;
end

vdaoengine.analysis.grammars.Utils.setNodesForGraph(mml.principalTree,NodePositions);
vdaoengine.analysis.grammars.Utils.setEdgesForGraph(mml.principalTree,e);

% nodenum = mml.principalTree.Nodes.size();
% for i=1:nodenum
%     display(mml.principalTree.Nodes.get(i).x);
% end

mml.principalTree.defineStarsFromPrimitiveGraphStructure();

mml.computeLayout();

NodePositions2D = mml.computedLayout.getNodePositions();
