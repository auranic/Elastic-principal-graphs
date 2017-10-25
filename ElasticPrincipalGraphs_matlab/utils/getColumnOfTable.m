function [strings] = getColumnOfTable(FileName,ColumnNumber)

Table = importdata(FileName);

for i=1:size(Table,1)
    classes = strsplit(char(Table(i)));
    class = classes(ColumnNumber);
    strings(i) = class;
end

end
