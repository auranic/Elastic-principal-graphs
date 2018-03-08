function [PenalizedEnergy] = ComputePenalizedElasticEnergy(data,graph,part)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

alpha = graph.BranchingControls(1);
beta = graph.BranchingControls(2);

    % Calculate MSE by usage dists
    MSE = part.dists' * data.Weights / sum(data.Weights);

    % Decompose ElasticMatrix
    Lambda = triu(graph.Lambdas,1);
    StarCenterIndices = find(graph.Mus>0);
    
    % Calculate edge potential
    [row,col] = find(Lambda);
    dev = graph.NodePositions(row, :)- graph.NodePositions(col, :);
    l = Lambda(Lambda>0);
    
    % Now lambdas are modified accordingly to the penalty
    Ks = sum(graph.Lambdas>0);
    lp = max(Ks(row),Ks(col));
    lp = lp-2;
    lp(lp<0) = 0; 
    lpenalized = l+alpha*lp';   %remove -1
    
    EP = sum(lpenalized(:) .* sum(dev.^2, 2));

    % Form indeces
    indL = (graph.Lambdas) > 0;

    % Calculate harmonicity potential
    RP = 0;
    for i=1:size(StarCenterIndices,1)
        leaves = indL(:,StarCenterIndices(i));
        K = sum(leaves);
        dev = graph.NodePositions(StarCenterIndices(i),:)...
            - sum(graph.NodePositions(leaves, :)) / K;
        %RP = RP + graph.Mus(StarCenterIndices(i)) * power(K,beta) * sum(dev .^ 2);
        CorrectedMu = graph.Mus(StarCenterIndices(i));
        if K>2 
            CorrectedMu = CorrectedMu/beta;
        end
        RP = RP + CorrectedMu  * sum(dev .^ 2);
    end
    % Total energy is the sum.
    PenalizedEnergy = MSE + EP + RP;



end

