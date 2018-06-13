clear all;

mat = load('tree_300.mat');

alpha = 0.005;
lambda = 0.03;

X = mat.X; 

Npoints = size(X,1)
dim = size(X,2)

maxNumNodes = 50;

figure;
[np,ed,ReportTable] = computeElasticPrincipalGraph(X,maxNumNodes,'BranchingControls',[alpha 1],'Lambda',lambda,'Plots',0); 
plot(X(:,1),X(:,2),'k.','MarkerSize',5); 
drawGraph2D(np,ed,'Color','r'); 
drawnow;

figure;
dX = X(randsample(Npoints,150),:);

[np,ed,ReportTable] = computeElasticPrincipalGraph(dX,maxNumNodes,'BranchingControls',[alpha 1],'Lambda',lambda,'Plots',0); 
plot(dX(:,1),dX(:,2),'k.','MarkerSize',5); 
drawGraph2D(np,ed,'Color','r'); 
drawnow;


    inflationFactor = 20;
    X1 = zeros(size(X,1)*inflationFactor,size(X,2));
    k=1;
    STDV = std(X);
    if 1
        for i=1:size(X,1)
            for j=1:inflationFactor
                r = rand(1,size(X,2));
                p = 0.4*STDV.*r;
                X1(k,:) = X(i,:)+p;
                k=k+1;
            end
        end
        X = X1;
    end

figure;
[np,ed,ReportTable] = computeElasticPrincipalGraph(X,maxNumNodes,'BranchingControls',[alpha 1],'Lambda',lambda,'Plots',0); 
plot(X(:,1),X(:,2),'k.','MarkerSize',5); 
drawGraph2D(np,ed,'Color','r'); 
drawnow;
