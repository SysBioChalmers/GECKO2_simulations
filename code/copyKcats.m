function newModel = copyKcats(model1,model2)
for i =1:length(model1.enzymes)
    enzyme = model1.enzymes{i};
    %Find enzyme in model2
    enzIdx  = find(contains(model2.metNames,enzyme),1);
    if ~isempty(enzIdx)
        enzRxns = find(model2.S(enzIdx,:));
        %Exclude enzyme usage reaction
        enzRxns = enzRxns(1:(end-1));
        %Get KCats from model2
        Kcats   = model2.S(enzIdx,enzRxns);
        %Now search enzyme in model1 
        enzIdx  = find(contains(model1.metNames,enzyme),1);
        enzRxns = find(model1.S(enzIdx,:));
        enzRxns = enzRxns(1:(end-1));
        model1.S(enzIdx,enzRxns) = Kcats;
    end 
end 
newModel = model1;
end