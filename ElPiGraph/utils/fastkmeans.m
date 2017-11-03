function centroids = fastkmeans(data, k, maxNumberOfIterations)
%This function present standard kmeans clusterisation with two possible
%methods of calculations.
%Input features
%   data is n-by-m matrix of data points. Each row perpresents one
%       observations. 
%   k is required number of centroids
%   meth is method of distances calculation
%   0 is standard usage of norm function
%   1 is usage of pdist2 function
%   2 is calculation by my function of distance calculation
%
%Return k-by-m array of centroids
%

    % Treshold of calculation splitting to not overflow memory in elemnts
    % (each element require 8 byte).
    MaxArray = 10000000; %10 mln. elements or 80M memory

    % Get data sizes
    n  = size(data, 1);
    
    % Calculate threshold as number of points
    maxPoints = floor(MaxArray / k);
    
    % To rovide better comparability I implements stupid initialisation by
    % the first k vectors
    centroids = data(1:k, :);
    
    % Vector of "current" associations. All points are not associated with
    % any centroids.
    ass = zeros(n, 1);
    
    %Main loop
    for i=1:maxNumberOfIterations
        % Copy old associations for comparison
        oldAss = ass;
        % Calculate new associations
        ass = associate(data, centroids, maxPoints);
        % Check stop condition
        %disp(n - sum(oldAss == ass));
        if ~any(oldAss - ass)
            % All associations are the same
            break;
        end;
        %Recalculte centroids I am not sure that it is fastest
        %implementation
        for j = 1:k
            centroids(j, :) = mean(data(ass == j, :));
        end
    end
end

function ass = associate(data, centroids, maxPoints)
% Calculate associations of data points from data with nearest centroids.
% We can calculate distances for maxPoints datapoints simultaneously.
% dataLength is vector of squared lengthes of data points. Moreover squared
% distance is ||x||^2+||y(i)||^2-2(x,y(i)) where ||x|| is length of data vector.
% Since we search min of distances from the same x to all centroids y(i)
% then we can omit summand ||x||^2.
    n = size(data, 1);
    ass = zeros(n, 1);
    %Calculate squared length of centroids
    cent = centroids';
    centrLength = sum(cent.^2);
    for i = 1:maxPoints:n
        % Define last element for calculation
        last = i + maxPoints - 1;
        if last > n
            last = n;
        end
        % Prepare index
        ind = i:last;
        % Calculate distances
        d = bsxfun(@minus, centrLength, 2 * (data(ind,:) * cent));
        [~, ass(ind)] = min(d,[],2);
    end
end