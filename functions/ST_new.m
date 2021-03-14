%force lactate 
model = readSBML('/Users/jessica/SystemsBiology/S_thermophillus_v1_03.xml',1000);
modelLB = readSBML('/Users/jessica/SystemsBiology/LBUL.xml',1000);
modelLB = addReaction(modelLB, 'EX_caspep_e', 'metaboliteList', {'caspep[e]'},...
'stoichCoeffList', [-1]);
model = removeRxns(model,{'PEPTabc','CASL','CASabc','PEPTL','EX_mpept_e_','EX_caspep_e_'});
model = addReactionFromModelRef(model,{'CASL'},modelLB);
model = addReactionFromModelRef(model,{'CASabc'},modelLB);
model = addReactionFromModelRef(model,{'EX_caspep_e'},modelLB);
%model.lb(587)=-1000;
model.lb(518)=1;

model.lb(14) =-1000; model.ub(14)= 1000; model.lb(312)=0.1; model.ub(312)=1000; %to force ST to produce fol - FOLt, EX_fol
media = getMediaFromModel(model);                   
nutrients = regexprep(media.reactions,{'EX_','_e'},{'',''});

fba = optimizeCbModel(model);
positions_exchange_rxns = find(startsWith(model.rxns,'EX_'));

pos_exchange_with_flux = intersect(find(fba.x~=0), positions_exchange_rxns);
exchange_rxns_with_flux = model.rxns(pos_exchange_with_flux);
products =  regexprep(setdiff(exchange_rxns_with_flux, media.reactions),{'EX_','_e'},{'',''});
metabolites_to_track = union(nutrients,products);
n_metabolites_to_track = length(metabolites_to_track);
exchangeRxns_metabolites_to_track = strcat('EX_', metabolites_to_track, '_e');
%exchangeRxns_metabolites_to_track(9,1)={'EX_caspep_e_'};
%model = removeRxns(model,rxnsToRemove_)
% iMM904d = addReactionFromModelRef(iMM904b,{'PTAr'},bigg);

n_rows = 5;
n_cols = 5;
grid0 = zeros(n_rows,n_cols);

initialConcentrations = cell(size(metabolites_to_track));
for i = 1:length(metabolites_to_track)
    initialConcentrations{i} = zeros(n_rows,n_rows);
end

initialbiomassConcentration = zeros(n_rows,n_rows);

% we put a little bit of nutrients and cells at position 50,50 at time 0.
row = 2;
col = 3;

for row = 1:1:3
    for col= 2:1:4
        
        for i = 1:length(metabolites_to_track)
            if ismember(metabolites_to_track(i),nutrients)
                initialConcentrations{i}(row,col) = 10; %what units?
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
metaboliteConcentrations = cell(n_metabolites_to_track,n_iterations);
for m = 1:n_metabolites_to_track
    for i = 2:n_iterations
        metaboliteConcentrations{m,i} = grid0;
    end
end

for i = 1:n_metabolites_to_track
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
            concentrations = zeros(n_metabolites_to_track,1);
            for m = 1:n_metabolites_to_track
                c = metaboliteConcentrations{m,n};
                concentrations(m) = c(i,j);
            end
            biomassConcentration = BiomassConcentrations{n}(i,j);
            disp('')
            if all(concentrations==0)
                
            else
                %model.lb(14) =-1000; model.ub(14)= 1000; model.lb(315)=0.1; model.ub(315)=1000;
                [concentrationMatrix, excRxnNames, timeVec,...
                    biomassVec] = dynamicFBA(model, exchangeRxns_metabolites_to_track,...
                    concentrations, biomassConcentration, dt, 1, exchangeRxns_metabolites_to_track);
                %ons =  getPosOfElementsInArray(exchangeRxns_metabolites_to_track, excRxnNames);
                 
                concentrationMatrix = full(concentrationMatrix); %25x
                BiomassConcentrations{n+1}(i,j) = biomassVec(end);
                for m = 1:n_metabolites_to_track
                    if ismember(exchangeRxns_metabolites_to_track(m),excRxnNames)
                        pos =  find(strcmp(excRxnNames,exchangeRxns_metabolites_to_track(m)));
                        metaboliteConcentrations{m,n+1}(i,j) = concentrationMatrix(pos,end);
                    end
                end
                
            end
            
            
        end
    end
end
