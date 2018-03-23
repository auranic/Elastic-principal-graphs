function drawPieChart(x, y, radius, Proportions, Colors)
%drawPieChart draw one pie chart with specified properties:
%   x and y are coordinates of centre of circle
%   radius is radius of circle
%   Proportions ia vector of couns of casses for each color
%   Colors is vector of colors

    % Trancform counts to proportions
    prop = Proportions / sum(Proportions);
    % Starts from zero angle
    angle = 0;
    % Draw chart
    if size(prop,2)==1
        % One class chart
        for i=1:size(Proportions,2)
            color = Colors(i,:);
            plot_circle(x,y,radius,color);
        end
    else
        % Several classes case
        for i=1:size(Proportions,2)
            th1 = angle;
            th2 = angle+prop(i)*2*pi;
            color = Colors(i,:);
            plot_arc(th1,th2,x,y,radius,color);
            angle = th2;
        end
    end
end


