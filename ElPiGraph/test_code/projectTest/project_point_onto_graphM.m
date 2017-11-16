function [MSE, X_projected, EdgeIndices, ProjectionValues] =...
    project_point_onto_graphM(X, node, Edges, partition)
%% This function calculates piece-wise linear projection of a dataset onto
%% the graph defined by NodePositions and Edges
% Input arguments:
%   X - is the n-by-m data matrix. Each row corresponds to one data point.
%   NodePositions - is k-by-m matrix of positions of the original graph
%       nodes in the same space as X.
%   Edges is k-by-2 matrix of integers. Edges(i,1) and Edges(i,2) specify
%       numbers of two vertex of i-th edge.
%   partition is n-by-1 vector. partition(i) is number of node which is
%       associated with data point X(i,:).
% 
% Outputs:
%   MSE is mean squared error (distance) from X to projection X onto graph
%   X_projected is matrix of the coordinates of data point projections
%   EdgeIndices is n-by-1 with number of edge which contains projection
%   ProjectionValues is n-by-1 real vector with values between 0 and 1
%       indicating where between Edge(index,1) and Edge(index,2) the
%       projection was done.  
%%
            %   Get edge list
            L1 = Edges(1, :);
            L2 = Edges(2, :);
            %1. get array of edged end
            V2 = node(L2,:);
            %2.	Form matrix of edge directions
            dir = node(L1, :) - V2;
            %3.	Calculate squared length of edge directions
            len = sqrt(sum(dir .^ 2, 2));
            % Normalize vector of directions
            dir = bsxfun(@rdivide, dir, len);
            % Calculate projections
            pr = bsxfun(@minus, X * dir', (V2 * dir')');
            % Copy projections to normalize
            prn = bsxfun(@rdivide, pr, len');
            % Non negativity
            prn(prn < 0) = 0;
            % Cut too long edges
            prn(prn > 1) = 1;
            %Final normalization of prn:
            prn = bsxfun(@times, prn, len');
            %Calculate distances:
            dist = bsxfun(@plus, sum(X .^ 2, 2), len')...
                + prn .^ 2 - 2 * X * V2' - 2 * prn .* pr;
            %Select the shortest distance:
            [dist, ind] = min(dist,[],2);
            % Form outputs
            MSE = mean(dist);
            X_projected = V2 + bsxfun(@times, ProjectionValues, dir(ProjectionValues, :)');
            EdgeIndices = ind;
            ProjectionValues = prn(:, ind);
end
