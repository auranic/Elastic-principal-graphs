function LabelColorMap = createLabelColorMapList(labels)
%createLabelColorMapList create default LabelColorMap on base of list of
%   labels.
%Input features
%   labels is vector of labels (one label for each observation).
%Return map of colours associated with labels

    % Create map
    labels = unique(labels);
    LabelColorMap = containers.Map;
    for i=1:length(labels)
        if i > 18
            LabelColorMap(char(labels(i))) = randomColor();
        else
            LabelColorMap(char(labels(i))) = predefinedColor(i);
        end
    end
end

function [color] = randomColor()
    %Create random color
    color(1) = rand();
    color(2) = rand();
    color(3) = rand();
end

function [color] = predefinedColor(i)
    sequence = [0.0, 0.0, 1.0;  0.0, 1.0, 0.0;  1.0, 0.0, 0.0;  
                0.0, 1.0, 1.0;  1.0, 0.0, 1.0;  0.8, 0.8, 0.0;  
                0.0, 0.0, 0.5;  0.0, 0.5, 0.0;  0.5, 0.0, 0.0;
                0.0, 0.5, 0.5;  0.5, 0.0, 0.5;  0.5, 0.5, 0.0;
                0.0, 0.5, 1.0;  0.0, 1.0, 0.5;  0.5, 0.0, 1.0;
                1.0, 0.0, 0.5;  0.5, 1.0, 0.0;  1.0, 0.5, 0.0];
    k = mod(i-1,size(sequence,1))+1;
    color = sequence(k,:);
end

