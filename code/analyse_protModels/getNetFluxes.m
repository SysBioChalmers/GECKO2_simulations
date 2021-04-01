function netFluxes = getNetFluxes(rxnIDs,fluxes,original_IDs)
current = pwd;
cd ..
%Get net fluxes for all reactions
netFluxes = zeros(length(original_IDs),1);
for i=1:length(original_IDs)
    ID = original_IDs{i};
    % Check if i-th reaction is reversible
    revFlag = any(find(contains(rxnIDs,ID) & contains(rxnIDs,'_REV')));
    % Map rxn indexes in flux dist table
    rxnIdxs = rxnMapping(ID,rxnIDs,revFlag);
    %calculate flux ratios and net flux in case the reaction is reversible
    if ~isempty(rxnIdxs)
        if length(rxnIdxs) == 2
            netFlux   = fluxes(rxnIdxs(1)) -fluxes(rxnIdxs(2));
        elseif length(rxnIdxs) == 1
            netFlux   = fluxes(rxnIdxs(1));
        else
            disp(ID)
            rxnIDs(rxnIdxs)
        end
        netFluxes(i) = netFlux;
    end
end
cd(current)
end