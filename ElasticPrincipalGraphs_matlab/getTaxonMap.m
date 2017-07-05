function [TaxonMap] = getTaxonMap(graph, data)
%Create TaxonMap to draw pie charts. Each element of map contains list of
%data points associated with this node
%Input features:
%   graph is object of vdaoengine.analysis.grammars.Graph class
%   data is n-by-m matrix data (each row of matrix contains one observation).
%   

    % Connection with Java to calculate
    javaclasspath({'VDAOEngine.jar'});
    % Get number of nodes
    numberOfNodes = graph.Nodes.size();
    % get sizes of data
    [numpoints, coordnum] = size(data);
    % Create object of vdaoengine.data.VDataSet class
    dataset = vdaoengine.data.VDataSet;
    dataset.massif = data;
    dataset.pointCount = numpoints;
    dataset.coordCount = coordnum;

    % Optimise graph
    elo = vdaoengine.analysis.grammars.ElasticEnergyOptimization(dataset, graph);
    elo.calcTaxons();
    % Create map
    TaxonMap = containers.Map;
    % Transform information to MatLab representation
    for i=1:numberOfNodes
        % Get number of data points which are associated with i-th node
        tx = elo.taxons.get(i-1);
        taxonSize = tx.size();
        % Create new array for data
        taxon = zeros(taxonSize,1);
        % Copy data from collection to array
        for j=1:taxonSize
            taxon(j) = tx.get(j-1)+1;
        end
        % Store new list to map
        TaxonMap(num2str(i)) = taxon;
    end
end