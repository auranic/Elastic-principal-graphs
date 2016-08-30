function [ output_args ] = drawPieChart(x,y,radius,Proportions,Colors)
%UNTITLED4 Summary of this function goes here
%   Proportions should be a row vector


prop = Proportions/sum(Proportions);

angle = 0;

if size(prop,2)==1 
   
for i=1:size(Proportions,2)
    color = Colors(i,:);
    plot_circle(x,y,radius,color); hold on;
end

    
else

for i=1:size(Proportions,2)
    th1 = angle;
    th2 = angle+prop(i)*2*pi;
    color = Colors(i,:);
    P = plot_arc(th1,th2,x,y,radius,color); hold on;
    angle = th2;
end

end

end


