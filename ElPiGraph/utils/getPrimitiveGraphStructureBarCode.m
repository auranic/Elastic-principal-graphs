function [barcode, N] = getPrimitiveGraphStructureBarCode(ElasticMatrix)
% getPrimitiveGraphStructureBarCode forms array N with number of stars with k
% leaves in element k and barcode.

    % Decompose ElasticMatrix
    L = ElasticMatrix - diag(diag(ElasticMatrix));
    Connectivities = sum(L>0);
    N = accumarray(Connectivities', 1);
    % Number of nodes
    barcode = ['|',int2str(size(ElasticMatrix,1))];
    % Continue by addition left part
    if length(N)<=2
        barcode = ['0|',barcode];
    else
        for i=3:length(N)
            barcode = [int2str(N(i)), '|', barcode]; %#ok<AGROW>
        end
    end
end

