function accuracyComplexityPlot(ReportTable, varargin)
%accuracyComplexityPlot draw plot of accuracy/complexity for selection of
%the optimal number of nodes
%ReportTable is output of ElPrincGraph or of some of higher level functions
%like computeElasticPrincipalCircle, computeElasticPrincipalCurve,
%computeElasticPrincipalGraph, or computeRobustElasticPrincipalGraph.
%
%Name,Value pairs allow specify one of accuracy measure and one of
%complexity measure:
%   'accuracy' can has following values
%       'MSE' is mean square error or assessment of data approximation
%           quality. 
%       'MSEP' is mean square error for projections to edge.
%       'FVE' (default) is fraction of explained variance. This value
%           always between 0 and 1. Greater value means higher quality of
%           data approximation. 
%       'FVEP' is fraction variance unexplained for projections to edge.
%   'complexity' can has following values
%       'UE' is elastic energy for edges stretching.
%       'UR' is elastic energy of deviation from harmonicity.
%       'URN' is UR * nodes 
%       'URN2' (default) is UR * nodes^2 

    complexity = 'URN2';
    accuracy = 'FVE';
    % Parse input arguments
    for i=1:2:length(varargin)
        if strcmpi(varargin{i},'accuracy')
            accuracy = varargin{i + 1};
        elseif strcmpi(varargin{i},'complexity')
            complexity = varargin{i + 1};
        end
    end
    

    % Extract data from structure
    try
        urn2 = ReportTable.(complexity);
    catch ME
        error(['Unknown complexity "', complexity, '".']);
    end
    try
        fvep = ReportTable.(accuracy);
    catch ME
        error(['Unknown accuracy "', accuracy, '".']);
    end
    if any(isnan(fvep))
        error(['Requested accuracy "', accuracy, '" was not really calculated.']);
    end
    bcodes = ReportTable.('BARCODE');
    % Create figure
    figure;
    % Draw graph
    plot(fvep, urn2, 'r-o', 'LineWidth', 2); 
    hold on;
    % Add markers of minimums
    for i = 2:size(urn2, 1) - 1
        xp = urn2(i - 1);
        x = urn2(i);
        xn = urn2(i + 1);
        if x < min(xp, xn)
            diff = abs(x - (xp + xn) / 2);
            if diff > 0.01
                s = bcodes(i,:);
                text(fvep(i), urn2(i), s, 'FontSize', 14);
            end
        end
    end
    title('Accuracy/Complexity plot', 'FontSize', 14);
    xlabel(['Accuracy: ', accuracy], 'FontSize', 12);
    ylabel(['Geometrical complexity:', complexity], 'FontSize', 12);
    rotate3d off;
end