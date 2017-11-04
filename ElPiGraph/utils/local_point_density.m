function [local_densities] = local_point_density(X, sampling)
%% calculates local density estimates for X in points defined by 
%% array of indices sampling

centers = X(sampling,:);
[partition,dists] = PartitionData(X,centers,100000,sum(X.^2,2));
npoints_in_centers = histc(partition,[min(partition):max(partition)]);

mean_dev_squared = mean(dists);

local_densities = zeros(length(sampling),1);

for i=1:length(sampling)
    if npoints_in_centers(i)>0
        local_densities(i) = sum(exp(-dists(partition==i)/mean_dev_squared));
    end
end

local_densities = local_densities/sum(local_densities);

end