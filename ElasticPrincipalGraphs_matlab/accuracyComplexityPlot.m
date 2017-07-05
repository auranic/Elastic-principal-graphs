function accuracyComplexityPlot(ReportTable)
%accuracyComplexityPlot draw plot of accuracy/complexity for selection of
%the optimal number of nodes
%ReportTable is output of computeElPT, computeElasticPrincipalTree, or
%   computeElasticPrincipalGraph.

    % Create figure
    figure;
    % Extract data from table
    urn2 = table2array(ReportTable(:,'URN2'));
    fvep = table2array(ReportTable(:,'FVE'));
    n = table2array(ReportTable(:,'NNODES'));
    % Draw graph
    plot(fvep,urn2,'r-o','LineWidth',2); 
    hold on;
    % Add markers of minimums
    for i=2:size(urn2,1)-1
        xp = urn2(i-1);
        x = urn2(i);
        xn = urn2(i+1);
        if x<min(xp,xn)
            diff = abs(x-(xp+xn)/2);
            if diff>0.01
                text(fvep(i),urn2(i),table2array(ReportTable(i,2)),'FontSize',14);
            end
        end
    end
    title('Accuracy/Complexity plot','FontSize',14);
    xlabel('Fraction of Explained variance','FontSize',12);
    ylabel('Geometrical complexity','FontSize',12);
    rotate3d off;
end