function [barcode] = getPrimitiveGraphStructureBarCode(ElasticMatrix)

Mu = diag(ElasticMatrix);
L = ElasticMatrix - diag(Mu);
Connectivities = sum(L>0);
maxc = max(Connectivities);
ee = 1:maxc;
N = histc(Connectivities,ee);
barcode = strcat('||',int2str(size(ElasticMatrix,1)));
if maxc<=2
    barcode = strcat('0',barcode);
else
    for i=3:size(N,2)
        if i~=size(N,2)
            barcode = strcat(strcat('|',int2str(N(i))),barcode);
        else
            barcode = strcat(int2str(N(i)),barcode);
        end
    end
end

end

