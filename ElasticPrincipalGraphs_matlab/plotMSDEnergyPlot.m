function plotMSDEnergyPlot(ReportTable,e)

dt = ReportTable.data;

n = dt(:,2);

energy = dt(:,1);
msep = dt(:,8);

h=plot(n,[msep,energy],'-','LineWidth',2); hold on;
legend(h,'MSD','El.Energy',-1);

epc1 = sum(e(2:end));
epc2 = sum(e(3:end));

plot([2 max(n)],[epc1 epc1],'k-');
plot([2 max(n)],[epc2 epc2],'m-');