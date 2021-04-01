%organism dependent variables
current = pwd;
orgCodes    = {'eco' 'yli' 'kma'};
model_src   = {'../../ecModels/iML1515.xml' ...
               '../Yarrowia_lipolytica_W29-GEM/ModelFiles/mat/iYali.mat' ...
               '../Kluyveromyces_marxianus-GEM/ModelFiles/mat/Kluyveromyces_marxianus-GEM.mat'};
ecModels = {'eciML1515' 'eciYali' 'eciSM966'};
for i=3:length(orgCodes)
    disp(['Running comparative FVA for: ' orgCodes{i} ])
    %load GEM
    if strcmpi(orgCodes(i),'eco')
        model = importModel(model_src{i});
    else
        load(model_src{i})
    end
    if isfield(model,'model')
        model = model.model;
    end
    %convert to RAVEN if necessary
    if isfield(model,'rules')
        model = ravenCobraWrapper(model);
    end
    %load ecModel
    load(['../ecModels/' ecModels{i} '/model/' ecModels{i} '_batch.mat'])
    cd(['../' orgCodes{i} '_scripts'])
    %Set media constraints
    parameters = getModelParameters;
    cSource    = parameters.c_source;
    if strcmpi(orgCodes(i),'kma')
       model = changeMedia_Original(model,[],[],1000);
    elseif strcmpi(orgCodes(i),'yli')
    	model = changeMedia_Original(model,'D-Glucose exchange','Min',1000);
    else
    	model = changeMedia_Original(model,'D-Glucose exchange');
    end
    ecModel_batch = changeMedia_batch(ecModel_batch,'D-Glucose exchange (reversible)');
    cd ../GECKO/geckomat/utilities/FVA
    [FVA_Dists_batch,indexes_batch,blocked_batch,stats_batch] = comparativeFVA(model,ecModel_batch,'D-Glucose exchange',false,1E-12,false);
    close all
    [FVA_Dists_chemo,indexes_chemo,blocked_chemo,stats_chemo] = comparativeFVA(model,ecModel_batch,'D-Glucose exchange',true,1E-12,false);
    close all
    cd(current)
    plotFVRcumDist(FVA_Dists_batch,'Batch conditions','FV range [mmol/gDw h]','relative frequency',[orgCodes{i} '_batch_FVA.tiff'])
    close all
    FVAtable = table(ecModel_batch.rxns(indexes_batch),ecModel_batch.rxnNames(indexes_batch),FVA_Dists_batch{1},FVA_Dists_batch{2});
    FVAtable.Properties.VariableNames = {'rxns' 'rxnNames' 'FV_GEM' 'FV_ecModel'};
    writetable(FVAtable,['../../results/' orgCodes{i} '_FVA_batch.txt'],'delimiter','\t','QuoteStrings',false)
    
    plotFVRcumDist(FVA_Dists_chemo,'chemostat conditions','FV range [mmol/gDw h]','relative frequency',[orgCodes{i} '_chemostat_FVA.tiff'])
    close all
    FVAtable = table(ecModel_batch.rxns(indexes_chemo),ecModel_batch.rxnNames(indexes_chemo),FVA_Dists_chemo{1},FVA_Dists_chemo{2});
    FVAtable.Properties.VariableNames = {'rxns' 'rxnNames' 'FV_GEM' 'FV_ecModel'};
    writetable(FVAtable,['../../results/' orgCodes{i} '_FVA_chemostat.txt'],'delimiter','\t','QuoteStrings',false)
    clc
end
