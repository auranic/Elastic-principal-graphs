%Test of kmeans

% Specify parameters
nObservations = 100000;
k = 10;
dim = 100;

% Generate data matrix
data = rand(nObservations, dim);

% Three startsa of kmeans
disp('Circle');
cent0 = kmeans(data, k, 0);
%cent0 = kmeansEps(data, k, 0);

disp('pdist2');
cent1 = kmeans(data, k, 1);
%cent1 = kmeansEps(data, k, 1);

disp('Dot product');
cent2 = kmeans(data, k, 2);
%cent2 = kmeansEps(data, k, 2);