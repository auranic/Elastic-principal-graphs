function accuracyComplexityPlot(ReportTable)

data = ReportTable.data;

urn2  = data(:,14);
fvep = data(:,9);
n = data(:,2);

plot(fvep,urn2,'r-o','LineWidth',2); hold on;


for i=2:size(urn2,1)-1
    xp = urn2(i-1);
    x = urn2(i);
    xn = urn2(i+1);
    if x<min(xp,xn)
        diff = abs(x-(xp+xn)/2);
        if diff>0.01
            %text(fvep(i),urn2(i),num2str(n(i)),'FontSize',14);
            text(fvep(i),urn2(i),ReportTable.textdata(i+1,2),'FontSize',14);
        end
    end
end