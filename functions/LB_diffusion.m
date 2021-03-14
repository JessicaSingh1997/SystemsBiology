modelLB = readSBML('/Users/jessica/SystemsBiology/LBUL.xml',1000);
mediaLB = getMediaFromModel(modelLB);                   
nutrientsLB = regexprep(mediaLB.reactions,{'EX_','_e'},{'',''});

fbaLB = optimizeCbModel(modelLB);

positions_exchange_rxnsLB = find(startsWith(modelLB.rxns,'EX_'));
exportSolutionToJson(modelLB,fbaLB.x,'test')
pos_exchange_with_fluxLB = intersect(find(fbaLB.x~=0), positions_exchange_rxnsLB);
exchange_rxns_with_fluxLB = modelLB.rxns(pos_exchange_with_fluxLB);
productsLB =  regexprep(setdiff(exchange_rxns_with_fluxLB, mediaLB.reactions),{'EX_','_e'},{'',''});
metabolites_to_trackLB = union(nutrientsLB,productsLB);
n_metabolites_to_trackLB = length(metabolites_to_trackLB);
exchangeRxns_metabolites_to_trackLB = strcat('EX_', metabolites_to_trackLB, '_e');
modelLB.lb(336)=0;
n_rows = 5;
n_cols = 5;
grid0 = zeros(n_rows,n_cols);

initialConcentrations = cell(size(metabolites_to_trackLB));
for i = 1:length(metabolites_to_trackLB)
    initialConcentrations{i} = zeros(n_rows,n_rows);
end

initialbiomassConcentration = zeros(n_rows,n_rows);

% we put a little bit of nutrients and cells at position 50,50 at time 0.
row = 2;
col = 3;

for row = 1:1:3
    for col= 2:1:4
        
        for i = 1:length(metabolites_to_trackLB)
            if ismember(metabolites_to_trackLB(i),nutrientsLB)
                initialConcentrations{i}(row,col) = 10; %what units?
                %initialConcentrations{18}(row,col)=100;
            else
                initialConcentrations{i}(row,col) = 0;
            end
        end
        
        initialbiomassConcentration(row,col) = 0.01; %what units?
    end
end



dt = 0.1; %hrs (user-defined parameter)
totalSimulationTime = 1; %hrs (user-defined parameter)
n_iterations = totalSimulationTime/dt;

%This is where the concentrations of metabolites over time and in space will be stored
% metaboliteConcentrations(i,j) is the grid0 at time point j, for metabolite
% i. metaboliteConcentrations(i,j) will be a matrix of size n_rows, n_cols.
metaboliteConcentrations = cell(n_metabolites_to_trackLB,n_iterations);
for m = 1:n_metabolites_to_trackLB
    for i = 2:n_iterations
        metaboliteConcentrations{m,i} = grid0;
    end
end

for i = 1:n_metabolites_to_trackLB
    metaboliteConcentrations{i,1} = initialConcentrations{i};
end
BiomassConcentrations = cell(1,n_iterations);
BiomassConcentrations{1} = initialbiomassConcentration;
for i = 2:n_iterations
    BiomassConcentrations{i} = grid0;
end

for n = 1:n_iterations-1
    for i = 1:n_rows
        for j = 1:n_cols
            concentrations = zeros(n_metabolites_to_trackLB,1);
            for m = 1:n_metabolites_to_trackLB
                c = metaboliteConcentrations{m,n};
                concentrations(m) = c(i,j);
            end
            biomassConcentration = BiomassConcentrations{n}(i,j);
            disp('')
            if all(concentrations==0)
                
            else
                %model.lb(14) =-1000; model.ub(14)= 1000; model.lb(315)=0.1; model.ub(315)=1000;
                [concentrationMatrix, excRxnNames, timeVec,...
                    biomassVec] = dynamicFBA(modelLB, exchangeRxns_metabolites_to_trackLB,...
                    concentrations, biomassConcentration, dt, 1, exchangeRxns_metabolites_to_trackLB);
                %ons =  getPosOfElementsInArray(exchangeRxns_metabolites_to_track, excRxnNames);
                 
                concentrationMatrix = full(concentrationMatrix); %25x
                BiomassConcentrations{n+1}(i,j) = biomassVec(end);
                for m = 1:n_metabolites_to_trackLB
                    if ismember(exchangeRxns_metabolites_to_trackLB(m),excRxnNames)
                        pos =  find(strcmp(excRxnNames,exchangeRxns_metabolites_to_trackLB(m)));
                        metaboliteConcentrations{m,n+1}(i,j) = concentrationMatrix(pos,end);
                    end
                end
                
            end
            
            
        end
    end
end

