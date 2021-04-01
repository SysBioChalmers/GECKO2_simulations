function [absUsages,relUsages]=getEnzUsageDist(model,fluxDist,enzymes)
if nargin<3
    enzymes   = model.enzymes;
end
absUsages = zeros(length(enzymes),1);
relUsages = ones(length(enzymes),1)*Inf;
for i=1:length(enzymes)
    enzRxn = find(contains(model.rxns,['prot_' enzymes{i}]));
    absUsages(i) = fluxDist(enzRxn);
    if model.ub(enzRxn)<1
        relUsages(i) = fluxDist(enzRxn)/model.ub(enzRxn);
    end
end
end