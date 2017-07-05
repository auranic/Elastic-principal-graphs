function [LabelColorMap] = createLabelColorMapList(labels)
%createLabelColorMapList create default LabelColorMap on base of list of
%   labels.
%Input features
%   labels is vector of labels (one label for each observation).
%Return map of colours associated with labels

    % Create map
    labels = unique(labels);
    LabelColorMap = containers.Map;
    for i=1:length(labels)
        LabelColorMap(char(labels(i))) = randomColor();
    end
end

function [color] = randomColor()
    %Create random color
    color(1) = rand();
    color(2) = rand();
    color(3) = rand();
end

