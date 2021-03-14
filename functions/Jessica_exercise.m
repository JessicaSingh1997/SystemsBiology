function Jessica_exercise

model = readSBML('/Users/jessica/SystemsBiology/iML1515.xml',1000);
media = getMediaFromModel(model);
nutrients = regexprep(media.reactions,{'EX_','_e'},{'',''});
model = changeRxnBounds(model,'EX_o2_e',0,'b');
nutrients = setdiff(nutrients,'o2');

fba = optimizeCbModel(model);
positions_exchange_rxns = find(startsWith(model.rxns,'EX_'));

pos_exchange_with_flux = intersect(find(fba.x~=0), positions_exchange_rxns);
exchange_rxns_with_flux = model.rxns(pos_exchange_with_flux);
products =  regexprep(setdiff(exchange_rxns_with_flux, media.reactions),{'EX_','_e'},{'',''});
metabolites_to_track = union(nutrients,products);
n_metabolites_to_track = length(metabolites_to_track);
exchangeRxns_metabolites_to_track = strcat('EX_', metabolites_to_track, '_e');

n_rows = 1;
n_cols = 2;
grid0 = zeros(n_rows,n_cols);

initialConcentrations = cell(size(metabolites_to_track));
for i = 1:length(metabolites_to_track)
    initialConcentrations{i} = zeros(n_rows,n_rows);
end

initialbiomassConcentration = zeros(n_rows,n_rows);

% we put a little bit of nutrients and cells at position 50,50 at time 0.
row = 50;
col = 50;

for row = 45:1:55
    for col= 45:1:55
        
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
totalSimulationTime = 5; %hrs (user-defined parameter)
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
                [concentrationMatrix, excRxnNames, timeVec,...
                    biomassVec] = dynamicFBA(model, exchangeRxns_metabolites_to_track,...
                    concentrations, biomassConcentration, dt, 1, exchangeRxns_metabolites_to_track);
                concentrationMatrix = full(concentrationMatrix);
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

% concentration of acetate (first position) at time 0 (first position), in
% location ((50,50)
m1 = metaboliteConcentrations{1,1}; %this is a matrix of 100x100
c1 = metaboliteConcentrations{1,1}(50,50);
% concentration of acetate (first position) at the end time point, in
% location ((50,50)
m2 = metaboliteConcentrations{1,end}; %this is a matrix of 100x100
%c2 = metaboliteConcentrations{1,end}(50,50);

map = [];
for i = 100:i
    map = [map; [0,0,i/100]];
    
end
colormap(map)
%acetate
figure
subplot(2,1,1)
c_t0 = metaboliteConcentrations{1,1};
imagesc(c_t0)
grid on
colormap(map)
subplot(2,1,2)
c_tend = metaboliteConcentrations{1,end};
imagesc(c_tend)
colormap(map)
grid on

%glucose
subplot(2,2,1)
c_t0 = metaboliteConcentrations{11,1};
subplot(2,2,2)
c_tend = metaboliteConcentrations{11,end};

%biomass
subplot(2,2,1)
c_t0 = BiomassConcentrations{1};
subplot(2,2,2)
c_tend = BiomassConcentrations{end};


end