function [TaxonMap] = getTaxonMap(graph, data)

javaclasspath({'VDAOEngine.jar'});

numberOfNodes = graph.Nodes.size();

		numpoints = size(data,1);
		coordnum = size(data,2);
		dataset = vdaoengine.data.VDataSet;
		dataset.massif = data;
		dataset.pointCount = numpoints;
		dataset.coordCount = coordnum;


elo = vdaoengine.analysis.grammars.ElasticEnergyOptimization(dataset, graph);

elo.calcTaxons();

TaxonMap = containers.Map;

for i=1:numberOfNodes
    
    tx = elo.taxons.get(i-1);
    taxonSize = tx.size();
    clear taxon;
    taxon = zeros(taxonSize,1);
    for j=1:taxonSize
        taxon(j) = tx.get(j-1)+1;
    end
    
    TaxonMap(num2str(i)) = taxon;
    
end