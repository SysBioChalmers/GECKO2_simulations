function plotCumDist(dataCell,titleStr,xStr,yStr,legendStrs,filterZeros)
if nargin<6
    filterZeros = false;
end
data     = [dataCell{1}, dataCell{2}, dataCell{3}];
nonZeros = find(sum(data,2)>1E-12);
%homogeneize legend strings
legendStrs = pad(legendStrs);
nonZeroDists = [];
figure
if length(dataCell)==3
    colors = {[0.1 0 0.8] [0.55 0.55 0.5] [0.8 0.6 0]};
end
for i=1:length(dataCell)
    dataX = dataCell{i};
    if filterZeros
        dataX = dataX(nonZeros);
    end
    nonZeroDists = [nonZeroDists,{dataX}];
    medianVal = median(dataX);
    legendStrs{i} = [legendStrs{i} ' (' num2str(numel(dataX)) ' / ' num2str(round(medianVal,4)) ')'];
    %perform statistical test with distribution #1
    if i>1
        [~,pval] = kstest2(nonZeroDists{1},nonZeroDists{i},'Alpha',0.01);
        if pval<=0.05
            legendStrs{i} = [legendStrs{i} ' *'];
            if pval<=0.01
                legendStrs{i} = [legendStrs{i} '*'];
            end
        end
    end
    [f, x] = ecdf(dataX);
    plot(x,f,'LineWidth',5,'Color',colors{i})
    set(gca, 'XScale', 'log','FontSize',18)
    xlim([0 1000])
    ylim([0 1])
    %axis(axisLimits)
    title(titleStr)
    xlabel(xStr)
    ylabel(yStr)
    hold on
end
legend(legendStrs,'location','southeast')
hold off
end