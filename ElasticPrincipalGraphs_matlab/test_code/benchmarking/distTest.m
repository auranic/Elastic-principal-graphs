%Test of kmeans

% Specify parameters
nObservations = 10000;
k = 10;
dim = 100;
N = 10;

% Generate data matrix
data = rand(nObservations, dim);

% Three startsa of kmeans
disp('Circle');
cent0 = associate(data, k, N, 0);

disp('pdist2');
cent1 = associate(data, k, N, 1);

disp('Dot product');
cent2 = associate(data, k, N, 2);
