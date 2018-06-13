NumberOfComputations = 100;

clear AllNodePositions;
clear AllAdjacencyMatrices;

%filename = 'C:\Datas\Peshkin\demo\day22\bootstrap\bs';
filename = 'C:\Datas\ElPiGraph\temp\tre23';
NumberOfNodes = 50;
SamplingPercentage = 0.1;
TrimmingRadius = 0.4;

N = size(data,1);
m = size(data,2);
K = ceil(N*(1-SamplingPercentage));

AllNodePositions = zeros(NumberOfNodes,m,NumberOfComputations);
AllAdjacencyMatrices = zeros(NumberOfNodes,NumberOfNodes,NumberOfComputations);


for i=1:NumberOfComputations

    display(sprintf('Iteration %i',i));
    
    sampling = randperm(N,K);
    factor = 1+(rand()-1)*0.4;
    tr = TrimmingRadius*factor;
    tic; [np,ed] = computeElasticPrincipalGraph(data(sampling,:),NumberOfNodes,'TrimmingRadius',tr,'Lambda',0.02,'Mu',0.1,'Plots',2); toc;
    drawnow;
    saveas(gcf,sprintf('%s%2i.png',filename,i));
    close all;
    
    AllNodePositions(:,:,i) = np(:,:);
    elm = MakeUniformElasticMatrix(ed,1,0);
    AllAdjacencyMatrices(:,:,i) = elm(:,:);

end

