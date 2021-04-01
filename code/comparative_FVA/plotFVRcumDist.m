function plotFVRcumDist(variable,titleStr,xStr,yStr,plotName)
figure
colors = {[0.1 0 0.8] [0.8 0.6 0]};
for i=1:length(variable)
    [f, x] = ecdf(variable{i});
    plot(x,f,'LineWidth',5,'Color',colors{i})
    set(gca, 'XScale', 'log','FontSize',24)
    xlim([1E-6 2E3])
    ylim([0 1])
    
    xticks([1E-6 1E-4 1E-2 1E0 1E2])
    %xticklabels({'-3\pi','-2\pi','-\pi','0','\pi','2\pi','3\pi'})
    yticks([0 0.25 0.5 0.75 1])
    grid on
    hold on
end
title(titleStr)
xlabel(xStr)
ylabel(yStr)
saveas(gcf,['../../results/Figure_3/' plotName])
hold off
end