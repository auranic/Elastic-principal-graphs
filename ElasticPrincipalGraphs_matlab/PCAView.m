function PCAView( Nodes, Edges, data, cpg, pc1, pc2, pc1FVE, pc2FVE )
%PCAView draw the dots from data and graph in the space of the two selected
%principal components. Data will be centralized before projection. Mean
%vector of data will be subtracted from nodes positions.
%
%Input features
%   Nodes is array of nodes positions
%   Edges is k-by-2 matrix of integers. Edges(i,1) and Edges(i,2) specify
%       numbers of two vertex of i-th edge.
%   data is n-by-m matrix data (each row of matrix contains one observation).
%   cpg is vdaoengine.analysis.grammars.ComputePrincipalGraph Java
%       container object cpg with various service functions. It is output
%       parameter of computeElPT, computeElasticPrincipalTree, or
%       computeElasticPrincipalGraph Description of this functions can be
%       found in ??? 
%   pc1 is optional argument: 
%       If it is omitted then 1 is used. 
%       If it is integer number then it is number of principal component to
%           use.
%       If it is m-by-1 array then it is vector to calculate projections
%   pc2 is optional argument: 
%       If it is omitted and pc1 is integer number then pc1 + 1 is used. 
%       If it is omitted and pc1 is vector then 1 is used. 
%       If it is integer number then it is number of principal component to
%           use.
%       If it is m-by-1 array then it is vector to calculate projections
%   pc1FVE is optional argument. Used as part of axis title if specified.
%       If pc1 is number of component then this value is calculated.
%   pc2FVE is optional argument. Used as part of axis title if specified.
%       If pc2 is number of component then this value is calculated.

    % Check input arguments
    if nargin < 4
        error('At least Nodes, Edges, data, and cpg must be specified');
    end
    % Centralise data
    means = mean(data);
    data = bsxfun(@minus, data, means);
    Nodes = bsxfun(@minus, Nodes, means);
    % Check pc1 and pc2
    if nargin < 5
        pc1 = 1;
        pc2 = 2;
        pc1FVE = 0;
        pc2FVE = 0;
    elseif nargin < 6
        if length(pc1) == 1
            pc2 = pc1 + 1;
            pc2FVE = 0;
        else
            pc2 = 1;
        end
    elseif nargin < 7
        pc1FVE = 0;
        pc2FVE = 0;
    elseif nargin < 8
        pc2FVE = 0;
    end
    
    % Calculate principal components if necessary
    if length(pc1) == 1 || length(pc2) == 1
        [vglobal, uglobal, explainedVariances] = pca(data_centered);
    end
    % Form coordinates to draw
    if length(pc1) == 1
        xData = uglobal(:,pc1);
        xNodes = Nodes*vglobal(:,pc1);
        pc1FVE = explainedVariances(pc1)/sum(explainedVariances);
    else
        xData = data*pc1(:);
        xNodes = Nodes*pc1(:);
    end
    if length(pc2) == 1
        yData = uglobal(:,pc2);
        yNodes = Nodes*vglobal(:,pc2);
        pc1FVE = explainedVariances(pc2)/sum(explainedVariances);
    else
        yData = data*pc2(:);
        yNodes = Nodes*pc2(:);
    end
    % Get size of nodes
    NodeSizes = cpg.graph.countNumberOfPointsProjected(cpg.dataset) + 1;
    % Create figure
    figure;
    % Draw tree with specified sizes of nodes
    drawGraph2D([xNodes,yNodes],Edges,'NodeSizes',NodeSizes);
    % Draw data
    plot(xData,yData,'ko','MarkerSize',2);
    title('PCA view of principal tree','FontSize',20);
    if pc1FVE>0
        xlabel(sprintf('PCx(%2.2f%%)',pc1FVE*100),'FontSize',20);
    else
        xlabel('PCx','FontSize',20);
    end
    if pc2FVE>0
        ylabel(sprintf('PCy(%2.2f%%)',pc2FVE*100),'FontSize',20);
    else
        ylabel('PCy','FontSize',20);
    end
end

