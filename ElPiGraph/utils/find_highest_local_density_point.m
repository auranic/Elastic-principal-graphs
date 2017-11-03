function [dense_point_index, neighbours] = find_highest_local_density_point(X, sampling)
%% finds a point in a dataset which probably marks the region with the highest density
%   

[local_densities] = local_point_density(X, sampling);

k = find(local_densities==max(local_densities));

dd = sum((X-repmat(X(sampling(k),:),size(X,1),1)).^2');

[dds,inds] = sort(dd);

kd = find(dds>0);

dense_point_index = sampling(k);
neighbours = inds(kd);

end

