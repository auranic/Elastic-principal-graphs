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
    centroids = data(1:k,:);
    
    
    % Vector of "current" associations. All points are not associated with
    % any centroids.
    ass = zeros(n,1);
    
    meth=2;
    
    if meth == 2
        % calculate squared length of all data points
        dataLength = sum(data.^2, 2);
    end
    
    %Main loop
    for i=1:maxNumberOfIterations
        % Copy old associations for comparison
        oldAss = ass;
        % Calculate new associations
        if meth == 0
            ass = associate0(data, centroids);
        elseif meth == 1
            ass = associate1(data, centroids, maxPoints);
        elseif meth == 2
            ass = associate2(data, centroids, maxPoints, dataLength);
        end
        % Check stop condition
        %disp(n - sum(oldAss == ass));
        if ~any(oldAss-ass)
            % All associations are the same
            break;
        end;
        %Recalculte centroids I am not sure that it is fastest
        %implementation
        for i = 1:k
            centroids(i,:) = mean(data(ass == i,:));
        end
    end

end

function ass = associate0(data, centroids)
    % Calculate associations of data points from data with nearest
    % centroids
    n = size(data, 1);
    k = size(centroids, 1);
    ass = zeros(n, 1);
    for i = 1:n
        % Loop for points
        best = Inf;
        bestInd = -1;
        for j = 1:k
            d = norm(data(i,:)-centroids(j,:));
            if d<best
                best = d;
                bestInd = j;
            end
        end
        ass(i) = bestInd;
    end
end

function ass = associate1(data, centroids, maxPoints)
    % Calculate associations of data points from data with nearest
    % centroids. We can calculate distances for maxPoints datapoints
    % simultaneously.
    n = size(data, 1);
    ass = zeros(n, 1);
    for i = 1:maxPoints:n
        % Define last element for calculation
        last = i + maxPoints - 1;
        if last > n
            last = n;
        end
        % Prepare index
        ind = i:last;
        % Calculate distances
        d = pdist2(data(ind,:),centroids);
        [~, ass(ind)] = min(d,[],2);
    end
end

function ass = associate2(data, centroids, maxPoints, dataLength)
    % Calculate associations of data points from data with nearest
    % centroids. We can calculate distances for maxPoints datapoints
    % simultaneously. dataLength is vector of squared lengthes of data
    % points. 
    
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
        d = bsxfun(@plus, dataLength(ind), centrLength) - 2 * (data(ind,:) * cent);
        [~, ass(ind)] = min(d,[],2);
    end
end