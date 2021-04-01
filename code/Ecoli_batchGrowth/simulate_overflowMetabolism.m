%simulate_overflowMetabolism(model)
git('clone https://github.com/SysBioChalmers/GECKO.git')
cd GECKO
git ('devel')
clc
cd ..
load('../model/eciML1515_batch.mat')
%Get relevant exchange reaction indexes
bio_idx  = find(contains(model.rxnNames,'biomass pseudoreaction'));
glc_idx  = find(contains(model.rxnNames,'D-glucose exchange (reversible)'));
EtOH_idx = find(strcmpi(model.rxnNames,'ethanol exchange'));
pyr_idx  = find(strcmpi(model.rxnNames,'pyruvate exchange'));
oxy_idx  = find(contains(model.rxnNames,'oxygen exchange (reversible)'));
CO2_idx  = find(strcmpi(model.rxnNames,'carbon dioxide exchange'));
ace_idx  = find(strcmpi(model.rxnNames,'acetate exchange'));
gly_idx  = find(strcmpi(model.rxnNames,'glycerol exchange'));
for_idx  = find(strcmpi(model.rxnNames,'formate exchange'));
indexes = [bio_idx glc_idx oxy_idx CO2_idx EtOH_idx ace_idx gly_idx for_idx pyr_idx];
%Get WT values
cd ../scripts
model = changeMedia_batch(ecModel_batch,'D-Glucose exchange (reversible)');
sol = solveLP(model,1);
WTgRate = sol.x(bio_idx);
%Initialize
exchanges = zeros(0,9);

cd ../studies/GECKO/geckomat/utilities
Drate = 0;
iterations = 100;
for i=1:100
    Drate = (WTgRate)*(i-1)/(iterations-1);
    temp = setParam(model,'lb',bio_idx,Drate);
    solution = simulateChemostat(model,Drate,[glc_idx bio_idx],true);
    printFluxes(temp,solution)
    exch_i = solution(indexes)';
    exchanges = [exchanges; exch_i];
end

plot(exchanges(:,1),exchanges(:,2),exchanges(:,1),exchanges(:,3),...
    exchanges(:,1),exchanges(:,4),exchanges(:,1),exchanges(:,5),...
    exchanges(:,1),exchanges(:,6),exchanges(:,1),exchanges(:,7),...
    exchanges(:,1),exchanges(:,8),exchanges(:,1),exchanges(:,9),'LineWidth',6)
legend({'Glucose' 'Oxygen' 'CO_{2}' 'Ethanol' 'Acetate' 'Glycerol' 'Formate' 'Pyruvate'})
    