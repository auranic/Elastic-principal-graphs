function plotMSDEnergyPlot(ReportTable,eigValues)
%plotMSDEnergyPlot draw graph of MSE and elastic energy optimization.
%ReportTable is output of computeElPT, computeElasticPrincipalTree, or
%   computeElasticPrincipalGraph.
%eigValues is optional attribute and represent array of eigenvalues of data
%   matrix.

    % Create figure
    figure;
    % Extract data from structure
    n = ReportTable.('NNODES');
    % Draw main graph
    h = plot(n,[ReportTable.('MSE'),...
        ReportTable.('ENERGY')],'-','LineWidth',2); 
    hold on;
    % Draw eigen values if presented
    if nargin == 2
        k = length(eigValues);
        epc1 = sum(eigValues(2:k));
        epc2 = sum(eigValues(3:k));
        plot([2 max(n)],[epc1 epc1],'k-');
        plot([2 max(n)],[epc2 epc2],'m-');
        legend('MSD','El.Energy',['Sum eigenvalues 2:',num2str(k)],...
            ['Sum eigenvalues 3:',num2str(k)]); 
    else
        legend(h,'MSD','El.Energy');
    end
    % Figure decoration
    title('MSD and Elastic energy','FontSize',14);
    xlabel('Number of nodes','FontSize',12);
    ylabel('MSD, Energy','FontSize',12);
    hold off;
end