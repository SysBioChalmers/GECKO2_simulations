%get_prot_models
%
% Script for generating condition-specific ecModels constrained with
% proteomics data for S. cerevisiae, Kluyveromyces marxianus, Yarrowia
% lipolytica grown at 0.1h-1 dilution rate on reference, high temperature,
% low pH and osmotic stress conditions.
%

for i=1:length(orgCode)
    orgCode= {'sce' 'kma' 'yli'};
    ecModel_names = {'ecYeastGEM' 'eciSM966' 'eciYali'};
    %number of experimental replicate per condition for each organism
    grouping = {[3 3 3 3] [3 3 3 3] [4 3 4]};
    %Replace orgnism-specific scripts in GECKO:
    scripts = dir([orgCode{i} '_scripts']);
    for j = 1:length(scripts)
        script = scripts(j).name;
        if contains(script,'.m') | contains(script,'.txt') | contains(script,'.mat') | contains(script,'.tsv')
            fullName   = [orgCode{i} '_scripts/' script];
            %Retrieve script path within GECKO
            GECKO_path = dir(['GECKO/**/' script]);
            if ~isempty(GECKO_path)
                GECKO_path = GECKO_path.folder;
                %Replace script in GECKO in its container subfolder
                copyfile(fullName,GECKO_path)
            end
        end
    end
    load(['ecModels/' ecModel_names{i} '/model/' ecModel_names{i} '_batch.mat'])
    load(['ecModels/' ecModel_names{i} '/model/' ecModel_names{i} '.mat'])
    ecModel = copyKcats(ecModel,ecModel_batch);
    cd GECKO/geckomat/utilities/integrate_proteomics
    generate_protModels(ecModel,grouping{i},[ecModel_names{i} '_prot'],ecModel_batch)
    cd ../../../..
    mkdir(['../ecModels/' ecModel_names{i} '/prot_constrained'])
    copyfile('GECKO/models/prot_constrained/*',['../ecModels/' ecModel_names{i} '/prot_constrained'])
    cd ../../../../
end