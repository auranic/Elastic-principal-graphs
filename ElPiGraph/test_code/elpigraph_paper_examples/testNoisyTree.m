data = load('tree23.data'); 
plot(data(:,1),data(:,2),'ko','MarkerSize',5);
drawnow;

robustness_radius = 0.4;
number_of_nodes = 30;
lambda = 0.02;
mu = 0.05;
alpha = 0.0;

display('Adding 20%% of noisy data points');

amount_of_noise = 20;
data_noise = add_uniform_background_noise(data,amount_of_noise);
[NodePositions,Edges] = computeRobustElasticPrincipalGraph(data_noise,number_of_nodes,robustness_radius,'Lambda',lambda,'Mu',mu,'BranchingControls',[alpha,1],'Plots',0);
plot(data_noise(:,1),data_noise(:,2),'k.','MarkerSize',5); hold on; 
plot(data(:,1),data(:,2),'r.','MarkerSize',5); 
drawGraph2D(NodePositions,Edges); 
title(sprintf('Percentage of noise points = %i%%',amount_of_noise));
drawnow;

display('Adding 100%% of noisy data points');

figure;

amount_of_noise = 100;
data_noise = add_uniform_background_noise(data,amount_of_noise);
[NodePositions,Edges] = computeRobustElasticPrincipalGraph(data_noise,number_of_nodes,robustness_radius*0.8,'Lambda',lambda,'Mu',mu,'BranchingControls',[alpha,1],'Plots',0);
plot(data_noise(:,1),data_noise(:,2),'k.','MarkerSize',5); hold on; 
plot(data(:,1),data(:,2),'r.','MarkerSize',5); 
drawGraph2D(NodePositions,Edges); 
title(sprintf('Percentage of noise points = %i%%',amount_of_noise));
drawnow;

display('Adding 1000%% of noisy data points');

figure;

amount_of_noise = 1000;
data_noise = add_uniform_background_noise(data,amount_of_noise);
[NodePositions,Edges] = computeRobustElasticPrincipalGraph(data_noise,number_of_nodes,robustness_radius*0.7,'Lambda',lambda,'Mu',mu,'BranchingControls',[alpha,1],'Plots',0);
plot(data_noise(:,1),data_noise(:,2),'k.','MarkerSize',5); hold on; 
plot(data(:,1),data(:,2),'r.','MarkerSize',5); 
drawGraph2D(NodePositions,Edges); 
title(sprintf('Percentage of noise points = %i%%',amount_of_noise));



