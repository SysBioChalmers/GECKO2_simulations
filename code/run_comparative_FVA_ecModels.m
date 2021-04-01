%organism dependent variables
orgCodes    = {'yeast' 'eco' 'yli' 'kma'};
model_src   = {'yeast-GEM/ModelFiles/mat/yeastGEM.mat' 'ecModels/iML1515.xml' ...
               'Yarrowia_lipolytica_W29-GEM/ModelFiles/mat/iYali.mat' ...
               'Kluyveromyces_marxianus-GEM/ModelFiles/mat/Kluyveromyces_marxianus-GEM.mat'};
ecModels = {'ecYeastGEM' 'eciML1515' 'eciYali' 'eciSM966'};
for i=1:length(orgCodes)
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
    parameters = getModelParameteres;
    cSource    = parameters.c_source;
    if strcmpi(orgCodes(i),'eco')
        idx = find(strcmp(model.rxnNames,'EX_glc__D_e'));
        model.rxnNames{idx} = 'D-Glucose exchange';
    end
    model   = changeMedia_Original(model,'D-Glucose exchange');
    ecModel_batch = changeMedia_batch(ecModel_batch,'D-Glucose exchange (reversible)');
    cd ../GECKO/geckomat/utilities/FVA
    [FVA_Dists_ecoB,indexes_ecoB,blocked_ecoB,stats_ecoB] = comparativeFVA(model,ecModel_batch,'D-Glucose exchange',false,1E-12,false);
    [FVA_Dists_ecoC,inde