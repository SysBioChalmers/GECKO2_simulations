%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Csources_simulations
%
% Cellular batch growth simulations on diverse carbon sources using
% ecYeastGEM
%
% Ivan Domenzain. Last modified: 2018-01-26
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('../../ecModels/ecYeastGEM/ecYeast7_GECKOv2.mat')
current = pwd;
figure
axis square
file_name = '../../data/growthRates_data_carbonSources.txt';
fID       = fopen(file_name);
data      = textscan(fID,'%s %s %f','delimiter','\t');
efe       = fclose('all'); 
media     = [{'YEP'}, {'MAA'}, {'Min'}];
count      = 1;
conditions = [];
gRateValues = [];
gRateValExps = [];
for i=1:length(media)
    media_indexes    = indexes_string(data{2},media{i},false);
    gRates_exp{1}{i} = data{1}(media_indexes);
    gRates_exp{2}{i} = data{3}(media_indexes);
    gRates_sim{i}    = [];
    SSres = 0;
    SStot = 0;
    for j=1:length(gRates_exp{2}{i})
        
        model = ecModel_batch;
        c_source = strcat(gRates_exp{1}{i}{j},' exchange (reversible)') 
        %c_source = strrep(c_source,'alpha,alpha-','');
        [model,pos]         = changeMedia_batch(model,c_source,media{i});
        gR_pos             = find(strcmpi(model.rxnNames,'growth'));
         model.c            = zeros(size(model.c));
         model.c(gR_pos)    = 1;
        
         solution           = solveLP(model);
         if solution.stat ~= 1
             disp(c_source)
             disp(solution.stat)
         end  
         model.lb(gR_pos)   = 0.999*solution.x(gR_pos);
         model.ub(gR_pos)   = solution.x(gR_pos);
         solution           = solveLP(model,1);
        gRates_sim{i}  = [gRates_sim{i};solution.x(gR_pos)];
        res                 = abs((gRates_exp{2}{i}(j)-gRates_sim{i}(j))/...
                              gRates_exp{2}{i}(j))*100;
        SSres               = SSres + res;
        

        switch gRates_exp{1}{i}{j}
                case 'D-fructose'
                    c_tag = 'Fru';
                case 'D-glucose'
                    c_tag = 'Glu';
                case 'sucrose'
                    c_tag = 'Suc';
                case 'maltose'
                   c_tag = 'Mal';
                case 'acetate'
                    c_tag = 'Ace';
                case 'D-galactose'
                    c_tag = 'Gal';
                case 'glycerol'
                    c_tag = 'Gly';
                case 'ethanol'
                    c_tag = 'Eth'; 
                case 'raffinose'
                    c_tag = 'Raf'; 
                case 'alpha,alpha-trehalose'
                    c_tag = 'Tre';
                case 'D-mannose'
                    c_tag = 'Man';

        end 
        flux_dist(:,count) = solution.x;
        count              = count+1;
        conditions         = [conditions;...
                              {char(strcat(media(i),string('_'),c_tag))}];
        gRateValues = [gRateValues;gRates_sim{i}(j)];
        gRateValExps = [gRateValExps;gRates_exp{2}{i}(j)];
        switch i
              case 1
                  marker = 'd';
              case 2
                  marker = 's';
              case 3
                  marker = 'o';
        end
        
        text(gRates_exp{2}{i}(j),gRates_sim{i}(j),c_tag,'FontSize',14)
        hold on

    end
    plot(gRates_exp{2}{i},gRates_sim{i},marker,'MarkerSize',15)
    hold on
    title('Max growth rate on different carbon sources','FontSize',30,'FontWeight','bold')
    ylabel('\mu_{max} predicted [h^{-1}]','FontSize',30,'FontWeight','bold');
    xlabel('\mu_{max} experimental [h^{-1}]','FontSize',30,'FontWeight','bold');
    RSq(i)       = (SSres)/length(gRates_exp{2}{i});
    legendStr(i) = strcat(media(i),' / e_{av}=',num2str(RSq(i)),'%');
end
x1 = linspace(0,1,1000);
plot(x1,x1)
hold on
legend(legendStr)
conditions = table(conditions,gRateValues,gRateValExps);
saveas(gcf,'../../results/Figure_1/ecYeastGEM_diverseCsources.jpg')
hold off
clc
disp('Resulting plot has been saved as: results/Figure_1/ecYeastGEM_diverseCsources.jpg')
clear
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function that receives a string and a cell array and returns the indexes
% in which the string appears on the array.
function matching = indexes_string(cell_array,str,flag)
    matching  = strfind(cell_array,str);
    if flag
       matching = find(~cellfun(@isempty,matching),1);
    else
        matching = find(~cellfun(@isempty,matching));
    end

end
