function [NodePositions,Edges,ReportTable,cpg,mml] = computeElasticPrincipalTree(data,NumNodes,varargin)

reduceDimension = 0;
newDimension = -1;

    for i=1:2:length(varargin)
        if strcmpi(varargin{i},'ReduceDimension')
            reduceDimension = 1;
            newDimension = varargin{i+1};
        elseif strcmpi(varargin{i},'')
            something = varargin{i+1};
        end
    end
    
mv = mean(data);    
mv1 = repmat(mv,size(data,1),1);
data_centered = data - mv1;

if reduceDimension
    [vglobal,uglobal,sglobal] = princomp(data_centered);
    if size(newDimension,2)>1
        startPC = newDimension(1);
        endPC = newDimension(2);
    else
        startPC = 1;
        endPC = newDimension(1);
    end
    perc = sum(sglobal(startPC:endPC))/sum(sglobal)*100;
    display(sprintf('Variance retained in %3.0f dimensions: %2.2f%%',(endPC-startPC+1),perc));
    data_centered = uglobal(:,startPC:endPC);
end

[NodePositions,Edges,ReportTable,cpg] = computeElPT(data_centered,NumNodes,varargin{:});

%%%%%%%%%%%%%%%%%%%%%%%%  Plots of MSE, elastic energy optimization
figure;
[v,u,explainedVariances] = princomp(data_centered);
if ~reduceDimension
    vglobal = v;
    uglobal = u;
end

plotMSDEnergyPlot(ReportTable,explainedVariances);
set(gcf,'Position',[11   338   612   219]);
title('MSD and Elastic energy','FontSize',14);
xlabel('Number of nodes','FontSize',12);
ylabel('MSD, Energy','FontSize',12);

%%%%%%%%%%%%%%%%%%%%%%%%  Accuracy/Complexity plot
figure;
accuracyComplexityPlot(ReportTable);
set(gcf,'Position',[11   665   612   316]);
title('Accuracy/Complexity plot','FontSize',14);
xlabel('Fraction of Explained variance','FontSize',12);
ylabel('Geometrical complexity','FontSize',12);

%%%%%%%%%%%%%%%%%%%%%%%% Show principal component view on principal tree and the data
figure;

NodeSizes = cpg.graph.countNumberOfPointsProjected(cpg.dataset);
ns = NodeSizes+1;

projectedNodesPCA = NodePositions*v;

if reduceDimension
varExpPC1 = sglobal(startPC)/sum(sglobal);
varExpPC2 = sglobal(startPC+1)/sum(sglobal);
else
varExpPC1 = explainedVariances(1)/sum(explainedVariances);
varExpPC2 = explainedVariances(2)/sum(explainedVariances);
end


plot(u(:,1),u(:,2),'ko','MarkerSize',2); hold on;
drawGraph2D(projectedNodesPCA(:,1:2),Edges,'NodeSizes',ns);
title('PCA view on principal tree','FontSize',20);
xlabel(sprintf('PC1(%2.2f%%)',varExpPC1*100),'FontSize',20);
ylabel(sprintf('PC2(%2.2f%%)',varExpPC2*100),'FontSize',20);
set(gcf,'Position',[641   339   695   643]);


%%%%%%%%%%%%%%%%%%%%%%%% Producing metro map layout
figure;

[NodePositions2D,mml] = computeMetroMapLayout(NodePositions,Edges); drawGraph2D(NodePositions2D,Edges,'NodeSizes',ns);
title('Metro map layout of the principal tree','FontSize',20);
set(gcf,'Position',[841   239   695   643]);

%%%%%%%%%%%%%%%%%%%%%%%% Preparing the output arguments

 if reduceDimension
     %Project nodes back into the initial, non-reduced space
     NodePositions = NodePositions*vglobal(:,startPC:endPC)';
 end
 
mv1 = repmat(mv,size(NodePositions,1),1);
NodePositions = NodePositions+mv1;