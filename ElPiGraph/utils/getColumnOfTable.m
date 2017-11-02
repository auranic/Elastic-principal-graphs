function [strings] = getColumnOfTable(FileName, ColumnNumber)
%getColumnOfTable load data from file to array. We assume that we read
%several rows of data with white space (space, tabulator, etc.) as
%separator. Cell array with content of the ColumnNumber column is returned.
    Table = importdata(FileName);
    N = size(Table,1);
    strings = cell(1,N);
    for i=1:N
        classes = strSplit(char(Table(i)));
        class = classes(ColumnNumber);
        strings(i) = class;
    end
end
