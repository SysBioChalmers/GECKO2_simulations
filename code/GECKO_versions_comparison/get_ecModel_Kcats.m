for dataset = {'old' 'new'}
load(['../../data/GECKO' dataset{1} '_enzData.mat'])
[m,n] = size(kcats.forw.kcats);
model = model_data.model;
values = [];
for i=1:n
    values = [values; kcats.forw.kcats(kcats.forw.kcats(:,i)>0,i)];
end
[m,n] = size(kcats.back.kcats);
for i=1:n
    values = [values; kcats.back.kcats(kcats.back.kcats(:,i)>0,i)];
end
T = table(values);
writetable(T,['../../results/' dataset{1} 'Kcats.txt'],'Delimiter','\t','QuoteStrings',false)

Kmatrix = kcats.tot.matrix;
vars = {'WC0' 'WC1' 'WC2' 'WC3'};
cmdStr = ['T_' dataset{1} '= table(Kmatrix(:,1),Kmatrix(:,2),Kmatrix(:,3),Kmatrix(:,4));'];
eval(cmdStr)
eval(['T_' dataset{1} '.Properties.VariableNames=vars;'])
end

WC_table = [];
for j=1:2
    Tvars = T_new;
    if j==2
        Tvars = T_old;
    end
    compTable = [];
    for i=1:width(T_new)
        compTable = [compTable; sum(table2array(Tvars(:,i)))];
    end
    WC_table = [WC_table,compTable];
end
WC_table = table(WC_table(:,1),WC_table(:,2),'VariableNames',{'GECKOv2' 'GECKO'});
writetable(WC_table,'../../results/WC_comparison.txt','Delimiter','\t','QuoteStrings',false)
writetable(T_new,'../../results/WC_GECKOv2.txt','Delimiter','\t','QuoteStrings',false)
clear