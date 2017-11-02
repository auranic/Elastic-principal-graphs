function [X] = ThreeFactorsMutualInhibitionAndCofactor(number_of_points,number_of_iterations, parameter_set)

tspan = [0;100];

% bifurcation 
if parameter_set==1
y0 = [5;1;1];
parameters = [1 0.1 5 1;1 0.1 5 1; 0.5 0.5 0 0];
noise_y0 = 0.1;
noise_parameters = 0.1;
end

X = [];

for i=1:number_of_iterations
y0i=y0.*(ones(size(y0,1),size(y0,2))+noise_y0*rand(size(y0,1),size(y0,2)));
parametersi = parameters.*(ones(size(parameters,1),size(parameters,2))+noise_parameters*rand(size(parameters,1),size(parameters,2)));
    
[TF_dynamics,t] = TF1_mutualinhibition_TF2_cofactorTF3(parametersi,tspan,y0i);

X = cat(1,X,TF_dynamics);

%h = plot(t,TF_dynamics,'o-'); legend(h,'TF1','TF2','TF3'); hold on;
%figure;
%plot3(TF_dynamics(:,1),TF_dynamics(:,2),TF_dynamics(:,3),'o-'); hold on;
%plot([min(TF_dynamics(:,1)) max(TF_dynamics(:,1))],[min(TF_dynamics(:,1)) max(TF_dynamics(:,1))],'b-'); hold on;
end

%downsampling
k = randperm(size(X,1));
X = X(k(1:number_of_points),:);

[v,u,s] = pca(X);
figure; 
%plot3(u(:,1),u(:,2),u(:,3),'o'); hold on;
plot(u(:,1),u(:,2),'o'); hold on;


end

