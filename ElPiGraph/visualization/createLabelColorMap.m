function [LabelColorMap] = createLabelColorMap(FileName, ColumnNumber)
% Create labelled colour map for labaels defined in column ColumnNumber in
% the text file FileName. All colours are randomly generated.

LabelColorMap = containers.Map;

labels = getColumnOfTable(FileName,ColumnNumber);

for i=1:size(labels,2)
    if ~LabelColorMap.isKey(labels(i))
        LabelColorMap(char(labels(i))) = randomColor();
    end
end

end

function [color] = randomColor()
    color(1) = rand();
    color(2) = rand();
    color(3) = rand();
end

