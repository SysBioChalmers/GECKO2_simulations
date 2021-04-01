function plotCumDist(dataCell,titleStr,xStr,yStr)
figure
for i=1:length(dataCell)
    dataX = dataCell{1};
    [f, x] = ecdf(dataX);
    plot(x,f,'LineWidth',5)
    set(gca, 'XScale', 'log')
    xlim([0 1000])
    ylim([0 1])
    %axis(axisLimits)
    title(titleStr)
    xlabel(xStr)
    ylabel(yStr)
    hold on
end
hold off
end