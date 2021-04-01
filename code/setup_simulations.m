%clone relevant repos
git clone https://github.com/SysBioChalmers/ecModels
cd ecModels
git('checkout 79805ff')
cd ..
git clone https://github.com/SysBioChalmers/GECKO
cd GECKO
git fetch --all --tags
git checkout tags/v2.0.1
cd ..

git('clone --depth 1 https://github.com/SysBioChalmers/Yarrowia_lipolytica_W29-GEM')
git('clone --depth 1 https://github.com/SysBioChalmers/yeast-GEM')
git('clone --depth 1 https://github.com/SysBioChalmers/Kluyveromyces_marxianus-GEM')
%organism dependent variables
orgCodes    = {'yeast' 'eco' 'yli' 'kma'};
model_src   = {'yeast-GEM/ModelFiles/mat/yeastGEM.mat' ...
               'Yarrowia_lipolytica_W29-GEM/ModelFiles/mat/iYali.mat' ...
               'Kluyveromyces_marxianus-GEM/ModelFiles/mat/Kluyveromyces_marxianus-GEM.mat'};

    
ecModel_src = {'ecModels/ecYeastGEM/model/ecYeastGEM_batch.mat' 'ecModels/eciML1515/model/eciML1515_batch.mat' ...
               'ecModels/eciYali/model/ecYeastGEM_batch.mat' ...
               'Kluyveromyces_marxianus-GEM/ModelFiles/mat/Kluyveromyces_marxianus-GEM.mat'};


