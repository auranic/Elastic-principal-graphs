function [data_noised] = add_uniform_background_noise(data,percentage)
%% Adding uniform noise to data background
%   
r = rand(floor(size(data,1)*percentage/100),size(data,2)); for i=1:size(data,2) r(:,i) = (max(data(:,i))-min(data(:,i)))*r(:,i)+min(data(:,i)); end;
data_noised = [data' r']';

end

