 %compare predictions
current = pwd;
orgs          = {'sce' 'yli' 'kma'};
ecModel_names = {'ecYeastGEM' 'eciYali' 'eciSM966'};
conditions    = {'Std' 'HiT' 'LpH' 'Osm'};
model_src     = {'yeast-GEM/ModelFiles/mat/yeastGEM.mat' ...
                 'Yarrowia_lipolytica_W29-GEM/ModelFiles/mat/iYali.mat' ...
                 'Kluyveromyces_marxianus-GEM/ModelFiles/mat/Kluyveromyces_marxianus-GEM.mat'};
exch_fluxes = [];
mkdir('../../results/Figure_3')
for i=1:length(orgs)
    parameters = [];
    disp(orgs{i})
    cd (current)
    %transfer org specific files to GECKo
    cd ..
    %Replace scripts in GECKO:
    disp(' ')
    disp('Transfering files to GECKO')
    fileNames = dir([orgs{i} '_scripts']);
    for j = 1:length(fileNames)
        fileName = fileNames(j).name;
        if ~startsWith(fileName,'.') & ~contains(fileName,'changeMedia_Original')
            fullName   = [[orgs{i} '_scripts/'] fileName];
            GECKO_path = dir(['GECKO/**/' fileName]);
            GECKO_path = GECKO_path.folder;
            copyfile(fullName,GECKO_path)
        end
    end
    %load gem
    load(model_src{i})
    if isfield('model','model')
        model = model.model;
    end
    %convert to RAVEN if necessary
    if isfield(model,'rules')
        model = ravenCobraWrapper(model);
    end
    %Load ecModel_batch
    load(['ecModels/' ecModel_names{i} '/model/' ecModel_names{i} '_batch.mat'])
    %load condition fermentation data
    fermData   = readtable([orgs{i} '_scripts/fermentationData.txt'],'delimiter','\t');
    byProducts = fermData.Properties.VariableNames(7:end);
    n          = height(fermData);
    %Initialize tables for saving results
    exch_fluxes_GEM = (cell(0,6));
    exch_fluxes_GEM = cell2table(exch_fluxes_GEM);
    exch_fluxes_GEM.Properties.VariableNames = {'condition' 'Drate' 'glucose' 'oxygen' 'CO2' 'NGAM'};
    exch_fluxes_ecM = exch_fluxes_GEM;
    exch_fluxes_pro = exch_fluxes_GEM;    
    subsystems = getSubSystem_str(model);
    fluxDist   = table(model.rxns,model.rxnNames,model.grRules,subsystems,'VariableNames',{'rxns' 'rxnNames' 'grRules' 'subSystems'});
    subsystems = mapEnzymeSubSystems(ecModel_batch.enzymes,ecModel_batch);
    absUsage   = table(ecModel_batch.enzymes,ecModel_batch.enzGenes,ecModel_batch.enzNames,subsystems,'VariableNames',{'enzymes' 'genes' 'enzNames' 'subSystems'});
    relUsage   = absUsage;
    significance_fluxDist = table({},{},{},{},'VariableNames',{'conditions' 'h' 'pval' 'model'});
    for k=1:n
        disp(conditions{k})  
        byP_fluxes = fermData{k,7:end};
        %Run a simulation with GEM
        cd([orgs{i} '_scripts'])
        c_source   = 'D-glucose exchange';
        glc_flux   = fermData{k,4};
        disp('*** GEM ')
        %set GUR and byprodcuts constraints according to data
        [tempGEM] = changeMedia_Original(model,c_source,'Min',glc_flux);
        for l=1:length(byProducts)
            index = find(strcmpi(tempGEM.rxnNames,[byProducts{l} ' exchange']));
            tempGEM.ub(index) = 1.05*byP_fluxes(l);
        end
        sol = solveLP(tempGEM,1);
        %Get predicted Drate
        objValue = min([-sol.f fermData.Drate(k)]);
        cd ../GECKO/geckomat
        parameters = getModelParameters;
        if k>1
            %fix dilution rate and set NGAM as objective to maximize
            bioIndx = find(model.c);
            tempGEM = setParam(tempGEM,'lb',bioIndx,0.9999*objValue);
            tempGEM = setParam(tempGEM,'ub',bioIndx,1.0001*objValue);
            tempGEM = setParam(tempGEM,'obj',parameters.NGAM,1);
            tempGEM = setParam(tempGEM,'lb',parameters.NGAM,0);
            tempGEM = setParam(tempGEM,'ub',parameters.NGAM,1000);
            sol1 = solveLP(tempGEM,1);
        else
            sol1 = sol;
        end
        if ~isempty(sol1.x)
            GEM_Drate = sol1.x(find(strcmpi(tempGEM.rxnNames,parameters.exch_names{1})));
            GEM_gluc  = sol1.x(find(strcmpi(tempGEM.rxnNames,'D-glucose exchange')));
            GEM_oxy   = sol1.x(find(strcmpi(tempGEM.rxnNames,'oxygen exchange')));
            GEM_CO2   = sol1.x(find(strcmpi(tempGEM.rxnNames,parameters.exch_names{4})));
            GEM_NGAM  = sol1.x(find(strcmpi(tempGEM.rxns,parameters.NGAM)));
            solution1 = sol1.x;
        else
            GEM_Drate = 0;
            GEM_oxy   = 0;
            GEM_CO2   = 0;
            GEM_NGAM  = 0;
            solution1 = zeros(length(tempGEM.rxns),1);
        end
        eval(['fluxDist.' conditions{k} '_GEM = solution1;']) 
        disp('*** ecModel ')
        cd limit_proteins
        P      = parameters.Ptot;%sumProtein(ecModel_batch);
        Pratio = fermData.Ptot(k)/fermData.Ptot(1);
        Ptot   = P*Pratio;
        if isfield('parameters','GAM')
            GAM = parameters.GAM;
        else 
            GAM = [];
        end
        %set exprimental GUR as constraints
        cd ../kcat_sensitivity_analysis
        if ~strcmpi(orgs{i},'sce')
            [temp_ecModel] = changeMedia_batch(ecModel_batch,'D-glucose exchange (reversible)','Min',fermData{k,4});
        else
            temp_ecModel = ecModel_batch;
        end
        for l=1:length(byProducts)
            index = find(strcmpi(temp_ecModel.rxnNames,[byProducts{l} ' exchange']));
            temp_ecModel.ub(index) = 1.05*byP_fluxes(l);
        end
        cSourceRxn   = find(strcmpi(temp_ecModel.rxnNames,parameters.exch_names{2}));
        temp_ecModel = setParam(temp_ecModel,'lb',cSourceRxn,0.99*fermData{k,4});
        temp_ecModel = setParam(temp_ecModel,'ub',cSourceRxn,1.01*fermData{k,4});
        %Get maximal Drate 
        sol      = solveLP(temp_ecModel);
        g_index  = find(temp_ecModel.c);
        objValue = min([-sol.f fermData.Drate(k)]);
        if k>1
            %Fix dilution rate
            temp_ecModel = setParam(temp_ecModel,'lb',g_index,0.99*objValue);
            temp_ecModel = setParam(temp_ecModel,'ub',g_index,1.01*objValue);
            %set NGAM as objective to maximize
            %temp_ecModel = setParam(temp_ecModel,'obj','prot_pool_exchange',-1);
            %temp_ecModel = setParam(temp_ecModel,'obj',parameters.NGAM,1);
            temp_ecModel = setParam(temp_ecModel,'ub',parameters.NGAM,1000);
            temp_ecModel = setParam(temp_ecModel,'lb',parameters.NGAM,0);
            sol2 = solveLP(temp_ecModel,1);
        else
            sol2 = sol;
        end
        if ~isempty(sol2.x)
            ecM_Drate = sol2.x(g_index);
            ecM_gluc  = sol2.x(find(strcmpi(temp_ecModel.rxnNames,parameters.exch_names{2})));
            ecM_oxy   = sol2.x(find(strcmpi(temp_ecModel.rxnNames,parameters.exch_names{3})));
            ecM_CO2   = sol2.x(find(strcmpi(temp_ecModel.rxnNames,parameters.exch_names{4})));
            ecM_NGAM  = sol2.x(find(strcmpi(temp_ecModel.rxns,parameters.NGAM)));
            solution2 = sol2.x;
        else
            ecM_Drate = 0;
            ecM_oxy   = 0;
            ecM_CO2   = 0;
            ecM_NGAM  = 0;
            solution2 = zeros(length(temp_ecModel.rxns),1);
        end
        %Store flux distribution
        cd (current)
        netFluxes = getNetFluxes(temp_ecModel.rxns,solution2,model.rxns);
        eval(['fluxDist.' conditions{k} '_ecM = netFluxes;']) 
        %Store enz usage distribution
        [absUsages,relUsages] = getEnzUsageDist(temp_ecModel,solution2,absUsage.enzymes);
        eval(['absUsage.' conditions{k} '_ecM = absUsages;']) 
        eval(['relUsage.' conditions{k} '_ecM = relUsages;'])
        
        disp(' ')
        disp('*** ecModel protConstrained')
        load(['../../ecModels/' ecModel_names{i} '_prot/' ecModel_names{i} '_prot_' conditions{k} '.mat'])
        [enzymes,iA,iB] = intersect(absUsage.enzymes,ecModelP.enzymes);
        absUsage  = absUsage(iA,:);
        relUsage  = relUsage(iA,:);
        
        pro_Drate = 0;
        pro_oxy   = 0;
        pro_CO2   = 0;
        pro_NGAM  = 0;
        sol3 = solveLP(ecModelP,1);
        if ~isempty(sol3.x)
            pro_Drate = sol3.x(find(strcmpi(ecModelP.rxnNames,parameters.exch_names{1})));
            pro_gluc  = sol3.x(find(strcmpi(ecModelP.rxnNames,parameters.exch_names{2})));
            pro_oxy   = sol3.x(find(strcmpi(ecModelP.rxnNames,parameters.exch_names{3})));
            pro_CO2   = sol3.x(find(strcmpi(ecModelP.rxnNames,parameters.exch_names{4})));
            pro_NGAM  = sol3.x(find(strcmpi(ecModelP.rxns,parameters.NGAM)));
            solution3 = sol3.x;
        end
        netFluxes = getNetFluxes(ecModelP.rxns,solution3,model.rxns);
        %Store flux distribution
        eval(['fluxDist.' conditions{k} '_ecP = netFluxes;']) 
        %Store enz usage distribution
        [absUsages,relUsages] = getEnzUsageDist(ecModelP,solution3,absUsage.enzymes);
        eval(['absUsage.' conditions{k} '_ecP = absUsages;']) 
        eval(['relUsage.' conditions{k} '_ecP = relUsages;'])   
        newRow1 = [conditions(k),{pro_Drate} {pro_gluc} {pro_oxy} {pro_CO2} {pro_NGAM}]; 
        exch_fluxes_pro = [exch_fluxes_pro; newRow1];
        newRow2 = [conditions(k),{ecM_Drate} {ecM_gluc} {ecM_oxy} {ecM_CO2} {ecM_NGAM}];
        exch_fluxes_ecM = [exch_fluxes_ecM; newRow2];
        newRow3 = [conditions(k),{GEM_Drate} {GEM_gluc} {GEM_oxy} {GEM_CO2} {GEM_NGAM}];
        exch_fluxes_GEM = [exch_fluxes_GEM; newRow3];
        clc
        %Obtain cumulative distributions for all flux dists
        fluxes_ecM = getNetFluxes(ecModel_batch.rxns,solution2,model.rxns);
        fluxes_ecP = getNetFluxes(ecModelP.rxns,solution3,model.rxns);
        cd ..
        plotCumDist({solution1 fluxes_ecM fluxes_ecP},[orgs{i} ' ' conditions{k}],'Flux [mmol/gDw h]','Cumulative frequency',{'GEM' 'ecM' 'ecMP'},true);
        saveas(gcf,['../results/Figure_S1/' orgs{i} '_' conditions{k} '_fluxDistsComp.tif'])
        close all
        %pairwiseComp({solution1 fluxes_ecM fluxes_ecP},{'GEM' 'ecM' 'ecMP'},[orgs{i} ' ' conditions{k}],'GEM Fluxes [mmol/gDw h]','ecModels fluxes [mmol/gDw h]',true);
        %close all
    end
    
    writetable(exch_fluxes_pro,['../results/' orgs{i} '_exchFluxes_pro.txt'],'delimiter','\t','QuoteStrings',false);    
    writetable(exch_fluxes_GEM,['../results/' orgs{i} '_exchFluxes_GEM.txt'],'delimiter','\t','QuoteStrings',false);
    writetable(exch_fluxes_ecM,['../results/' orgs{i} '_exchFluxes_ecM.txt'],'delimiter','\t','QuoteStrings',false);
    
    writetable(fluxDist,['../results/' orgs{i} '_fluxDist.txt'],'delimiter','\t','QuoteStrings',false);    
    writetable(absUsage,['../results/' orgs{i} '_absUsage.txt'],'delimiter','\t','QuoteStrings',false); 
    writetable(relUsage,['../results/' orgs{i} '_relUsage.txt'],'delimiter','\t','QuoteStrings',false); 
    disp(' ')
    close all
end           

