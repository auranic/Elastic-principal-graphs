data = load('test_code/branchingTest/thick_turn.data');

%tab = importdata('test_data\LLE_Pinello\Data_Jason.txt');
%data = tab.data;

nnodes = 50; computeElasticPrincipalGraph(data,nnodes,'Plots',2); 
x = get(gcf,'Position'); 
set(gcf,'Position',[x(1)-x(3) x(2) x(3) x(4)]); 
computeElasticPrincipalGraph(data,nnodes,'BranchingControls',[0.01,0],'Plots',2);

return;

%tab = importdata('test_data\LLE_Pinello\Data_Gottgens.txt');
tab = importdata('test_data\LLE_Pinello\Data_Jason.txt');

data = tab.data;

%[n,e] = computeElasticPrincipalGraph(data,30,'Plots',2,'Lambda',0.05,'Mu',0.1,'TrimmingRadius',0.02);
[n,e] = computeElasticPrincipalGraph(data,50,'Plots',2,'Lambda',0.01,'Mu',0.1);

figure;

plot3(n(:,1),n(:,2),n(:,3),'ko'); hold on; plot3(data(:,1),data(:,2),data(:,3),'r.');
for i=1:size(e,1)
    plot3([n(e(i,1),1) n(e(i,2),1)],[n(e(i,1),2) n(e(i,2),2)],[n(e(i,1),3) n(e(i,2),3)],'b-','LineWidth',3);
end