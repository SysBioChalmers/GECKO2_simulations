function [model,pos] = changeMedia_Original(model,c_source,media,flux)
% changeMedia_batch
%
% function that modifies aGEM and makes it suitable for batch growth
% simulations on different carbon sources.
%
% model:  An enzyme constrained model
% meadia: Media type ('YEP' for complex, 'MAA' minimal with Aminoacids,
%                          'Min' for minimal media)
% flux:   (Optional) A cell array with measured uptake fluxes in mmol/gDwh
%
% model: a constrained ecModel
%
% usage: [model,pos] = changeMedia_batch(model,c_source,media,flux)
%
% Ivan Domenzain        2019-02-13

% Give the carbon source (c_source) input variable with the following
% format: c_source  = 'D-glucose exchange'

if nargin<4
    flux = 1;
    if nargin<3
        media ='Min';
    end
end
%First allow any secretion and block all uptakes
[~,exchange]       = getExchangeRxns(model);
model.ub(exchange) = +1000;
model.lb(exchange) = 0;
pos                = getComponentIndexes(model,c_source);
%Block O2 and glucose production (avoids multiple solutions):
model.ub(strcmp(model.rxnNames,'oxygen exchange'))    = 0;
model.ub(strcmp(model.rxnNames,'D-glucose exchange')) = 0;

%The media will define which rxns to fix:
if strcmpi(media,'YEP')
    N = 25;     %Aminoacids + Nucleotides
elseif strcmpi(media,'MAA')
    N = 21;     %Aminoacids
elseif strcmpi(media,'Min')
    N = 1;      %Only the carbon source
end
%UB parameter (manually optimized for glucose on Min+AA):
b = 0.08;
%UB parameter (manually optimized for glucose complex media):
c = 2;
%Define fluxes in case of ec model:
if nargin < 4   %Limited protein    
    if N>1
       flux    = -b*ones(1,N);
       if N>21
           flux(22:25) = -c;
       end
    end
    flux(1) = -1000;
end
%Fix values as UBs:
for i = 1:N
    model.lb(pos(i)) = -flux(i);
end
%Allow biomass production
model.ub(find(model.c)) = 1000;
%Allow uptake of essential components
model = setParam(model, 'lb', 'y001654', -1000); % 'ammonium exchange';
model = setParam(model, 'lb', 'y002100', -1000); % 'water exchange' ;
model = setParam(model, 'lb', 'y001861', -1000); % 'iron(2+) exchange';
model = setParam(model, 'lb', 'y001992', -1000); % 'oxygen exchange';
model = setParam(model, 'lb', 'y002005', -1000); % 'phosphate exchange';
model = setParam(model, 'lb', 'y002060', -1000); % 'sulphate exchange';
model = setParam(model, 'lb', 'y001832', -1000); % 'H+ exchange' ;
model = setParam(model, 'lb', 'y001671', -1000); % Biotin . exchange
model = setParam(model, 'lb', 'y001548', -1000); % pantothenate exchange
model = setParam(model, 'lb', 'y001967', -1000); % Niconitate exchange
model = setParam(model, 'lb', 'y001947', -1000); % Myo-inositol
model = setParam(model, 'lb', 'y002067', -1000); % Thiamin (1+) exchange
model = setParam(model, 'lb', 'y002028', -1000); % Pyridoxine exchange
model = setParam(model, 'lb', 'y001604', -1000); % Aminobenzoic acid
%Block bicarbonate uptake
model = setParam(model, 'ub', 'y001663', 0); % 'bicarbonate production' ;


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pos = getComponentIndexes(model,c_source)
    pos(1)  = find(strcmpi(model.rxnNames,c_source))
    pos(2)  = find(strcmpi(model.rxnNames,'L-alanine exchange'));
    pos(3)  = find(strcmpi(model.rxnNames,'L-arginine exchange'));
    pos(4)  = find(strcmpi(model.rxnNames,'L-asparagine exchange'));
    pos(5)  = find(strcmpi(model.rxnNames,'L-aspartate exchange'));
    pos(6)  = find(strcmpi(model.rxnNames,'L-cysteine exchange'));
    pos(7)  = find(strcmpi(model.rxnNames,'L-glutamine exchange'));
    pos(8)  = find(strcmpi(model.rxnNames,'L-glutamate exchange'));
    pos(9)  = find(strcmpi(model.rxnNames,'glycine exchange'));
    pos(10) = find(strcmpi(model.rxnNames,'L-histidine exchange'));
    pos(11) = find(strcmpi(model.rxnNames,'L-isoleucine exchange'));
    pos(12) = find(strcmpi(model.rxnNames,'L-leucine exchange'));
    pos(13) = find(strcmpi(model.rxnNames,'L-lysine exchange'));
    pos(14) = find(strcmpi(model.rxnNames,'L-methionine exchange'));
    pos(15) = find(strcmpi(model.rxnNames,'L-phenylalanine exchange'));
    pos(16) = find(strcmpi(model.rxnNames,'L-proline exchange'));
    pos(17) = find(strcmpi(model.rxnNames,'L-serine exchange'));
    pos(18) = find(strcmpi(model.rxnNames,'L-threonine exchange'));
    pos(19) = find(strcmpi(model.rxnNames,'L-tryptophan exchange'));
    pos(20) = find(strcmpi(model.rxnNames,'L-tyrosine exchange'));
    pos(21) = find(strcmpi(model.rxnNames,'L-valine exchange'));
    pos(22) = find(strcmpi(model.rxnNames,'D-glucose exchange'));
end
