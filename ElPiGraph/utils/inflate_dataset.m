function [X_inflated] = inflate_dataset(X,inflationFactor,thickness)

X1 = zeros(size(X,1)*inflationFactor,size(X,2));
k=1;
STDV = std(X);

for i=1:size(X,1)
    for j=1:inflationFactor
        r = rand(1,size(X,2));
        p = thickness*STDV.*r;
        X1(k,:) = X(i,:)+p;
        k=k+1;
    end
end
X_inflated = X1;



end