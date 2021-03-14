%%%

%Set up Grid Space - GridR/L-final concentrations GridTempR/L where changes are tracked

gridLB = zeros(15,15); %x and y coordinate of the surface where the LB is - Right
gridTempLB=zeros(15,15); %where changes in concentration are stored for LB (Right)
gridST = zeros(15,15); %x and y coordinate of the surface where the ST is - Left
gridTempST=zeros(15,15); %where changes in concentration are stored for ST (Left)


%Initialize the position of Biomass on Grid
%Setting distance between cells

[j,k]= size(gridLB);
%Postion of LB biomass at t=0
idxR=ceil(j/2);
idyR=8;%ceil(3*k/4);
%Position of ST biomass at t=0
idxL=ceil(j/2);
idyL=7;%ceil(k/4);
gridDBR = gridLB;
gridDBL = gridST;
gridST(idxL,idyL)=0; %Initialize concentration of CF metabolite - Folic Acid to zero
gridLB(idxR,idyR)=0; %Initialize concentration of CF metabolite - mpept to zero


%Number of timesteps 
dt = 0.1; %s 
totalSimulationTime = 2; %s
n_iterations = totalSimulationTime/dt;


%%% TOY MODEL ST %%%

load('tmST.mat')
fbaST = optimizeCbModel(model);
mediaST = [{'EX_a_e'},{'EX_mpept_e'}];
nutrientsST = regexprep(mediaST,{'EX_','_e'},{'',''});
positions_exchange_rxnsST = find(startsWith(model.rxns,'EX_'));
pos_exchange_with_fluxST = intersect(find(fbaST.x~=0), positions_exchange_rxnsST);
exchange_rxns_with_fluxST = model.rxns(pos_exchange_with_fluxST);
productsST =  [{'c'},  {'d'}];
metabolites_to_trackST = union(nutrientsST,productsST);
n_metabolites_to_trackST = length(metabolites_to_trackST);
exchangeRxns_metabolites_to_trackST = strcat('EX_', metabolites_to_trackST, '_e');
n_rows = 1;
n_cols = 1;
grid0ST = zeros(n_rows,n_cols);

initialConcentrationsST = cell(size(metabolites_to_trackST));
for i = 1:length(metabolites_to_trackST)
    initialConcentrationsST{i} = 0;
end

initialbiomassConcentrationST= 0.01;
metaboliteConcentrationsST = cell(n_metabolites_to_trackST,n_iterations,1);
for m = 1:n_metabolites_to_trackST
    for i = 2:n_iterations
        metaboliteConcentrationsST{m,i}= grid0ST;
    end
end

BiomassConcentrationsST = cell(1,n_iterations);
BiomassConcentrationsST{1} = initialbiomassConcentrationST;
for i = 2:n_iterations
    BiomassConcentrationsST{i} = grid0ST;
end

%%% TOY MODEL LB %%%
load('tmLB')
fbaLB = optimizeCbModel(modelLB);
mediaLB = [{'EX_a_e'    },{'EX_f_e'}];
nutrientsLB = regexprep(mediaLB,{'EX_','_e'},{'',''});
positions_exchange_rxnsLB = find(startsWith(modelLB.rxns,'EX_'));
pos_exchange_with_fluxLB = intersect(find(fbaLB.x~=0), positions_exchange_rxnsLB);
exchange_rxns_with_fluxLB = modelLB.rxns(pos_exchange_with_fluxLB);
productsLB=[ {'mpept'    },  {'e'    }];
metabolites_to_trackLB = [{'mpept'};{'f'};{'a'};{'e'}];
n_metabolites_to_trackLB = length(metabolites_to_trackLB);
exchangeRxns_metabolites_to_trackLB = strcat('EX_', metabolites_to_trackLB, '_e');
grid0LB = zeros(n_rows,n_cols);

initialConcentrationsLB = cell(size(metabolites_to_trackLB));
for i = 1:length(metabolites_to_trackLB)
    initialConcentrationsLB{i} = zeros(n_rows,n_cols);
end
        
initialbiomassConcentrationLB= 0.01; %what units?
metaboliteConcentrationsLB = cell(n_metabolites_to_trackLB,n_iterations,1);
for m = 1:n_metabolites_to_trackLB
    for i = 2:n_iterations
        metaboliteConcentrationsLB{m,i}= grid0LB;
    end
end

BiomassConcentrationsLB = cell(1,n_iterations);
BiomassConcentrationsLB{1} = initialbiomassConcentrationLB;
for i = 2:n_iterations
    BiomassConcentrationsLB{i} = grid0LB;
end

%Change ST media concentrations

for p = 1:length(metabolites_to_trackST)
    if ismember(metabolites_to_trackST(p),nutrientsST)
            initialConcentrationsST{1} = 0.01;
            initialConcentrationsST{4} = 0.001;
    else
        initialConcentrationsST{p}= 0;
    end
end
for p = 1:n_metabolites_to_trackST
metaboliteConcentrationsST{p,1} = initialConcentrationsST{p};
end

%Change LB media concentrations

for q = 1:length(metabolites_to_trackLB)
    if ismember(metabolites_to_trackLB(q),nutrientsLB)
        initialConcentrationsLB{2}= 0;
        initialConcentrationsLB{3}= 0.01 ;%what units?
    else
        initialConcentrationsLB{q} = 0;
    end
end
for q = 1:n_metabolites_to_trackLB
    metaboliteConcentrationsLB{q,1} = initialConcentrationsLB{q};
end


%Diffusion component iterations 
diffusedMPEPT=0;
diffusedAccumulatedMPEPT=0;
diffusedFA=0;
diffusedAccumulatedFA=0;
folicAcidCF=0;
mpeptCF=0;

for i = 1:n_iterations
    diffusedAccumulatedMPEPT=0;
    diffusedAccumulatedFA=0;
    gridTempLB=gridLB;
    gridTempST=gridST;
    
        %Change LB media concentrations

    for q = 1:length(metabolites_to_trackLB)
        if ismember(metabolites_to_trackLB(q),nutrientsLB)
            initialConcentrationsLB{2}= folicAcidCF;
            initialConcentrationsLB{3}= 0.01 ;%what units?
        else
            initialConcentrationsLB{q} = 0;
        end
    end
    for q = 1:n_metabolites_to_trackLB
        if ismember(metabolites_to_trackLB(q),nutrientsLB)
            metaboliteConcentrationsLB{2,i} = initialConcentrationsLB{2};
        else
        metaboliteConcentrationsLB{q,1} = initialConcentrationsLB{q};
        end
    end

    if i~=n_iterations
    [folicacidProduced,metaboliteConcentrationsST,BiomassConcentrationsST] = dynamicST(i,n_rows, n_cols, n_metabolites_to_trackST,metaboliteConcentrationsST,BiomassConcentrationsST,model,dt,exchangeRxns_metabolites_to_trackST,mpeptCF);
    %[mpeptProduced,metaboliteConcentrationsLB,BiomassConcentrationsLB] = dynamicLB(i,n_rows, n_cols, n_metabolites_to_trackLB,metaboliteConcentrationsLB,BiomassConcentrationsLB, modelLB,dt,exchangeRxns_metabolites_to_trackLB,folicAcidCF);
    [mpeptProduced,metaboliteConcentrationsLB,BiomassConcentrationsLB] = dynamicLB(i,n_rows, n_cols, n_metabolites_to_trackLB,metaboliteConcentrationsLB,BiomassConcentrationsLB, modelLB,dt,exchangeRxns_metabolites_to_trackLB);
    end

    gridST(idxL,idyL)=gridST(idxL,idyL)+folicacidProduced; %ST - produces Folic Acid 
    gridLB(idxR,idyR)=gridLB(idxR,idyR)+mpeptProduced;%LB - produces mpept

    %MPEPT diffusion 
    for ix = 1:j
        for iy=1:k
            diffusedAccumulatedMPEPT=0;
            if(gridTempLB(ix,iy)==0)
                continue;
            end
            if(isIndexValid(gridLB,ix,iy-1))
                if(gridTempLB(ix,iy-1)<gridTempLB(ix,iy))
                    diffusedMPEPT=diffusion1(gridTempLB(ix,iy-1),gridTempLB(ix,iy),0.01,0.1);
                    gridLB(ix,iy-1)=gridLB(ix,iy-1)+diffusedMPEPT;
                    diffusedAccumulatedMPEPT=diffusedAccumulatedMPEPT+diffusedMPEPT;
                end
            end

            if(isIndexValid(gridLB,ix,iy+1))
                if(gridTempLB(ix,iy+1)<gridTempLB(ix,iy))
                    diffusedMPEPT=diffusion1(gridTempLB(ix,iy+1),gridTempLB(ix,iy),0.01,0.1);
                    gridLB(ix,iy+1)=gridLB(ix,iy+1)+diffusedMPEPT;
                    diffusedAccumulatedMPEPT=diffusedAccumulatedMPEPT+diffusedMPEPT;
                end
            end

            if(isIndexValid(gridLB,ix-1,iy))
                if(gridTempLB(ix-1,iy)<gridTempLB(ix,iy))
                    diffusedMPEPT=diffusion1(gridTempLB(ix-1,iy),gridTempLB(ix,iy),0.01,0.1);
                    gridLB(ix-1,iy)=gridLB(ix-1,iy)+diffusedMPEPT;
                    diffusedAccumulatedMPEPT=diffusedAccumulatedMPEPT+diffusedMPEPT;
                end
            end

            if(isIndexValid(gridLB,ix+1,iy))
                if(gridTempLB(ix+1,iy)<gridTempLB(ix,iy))
                    diffusedMPEPT=diffusion1(gridTempLB(ix+1,iy),gridTempLB(ix,iy),0.01,0.1);
                    gridLB(ix+1,iy)=gridLB(ix+1,iy)+diffusedMPEPT;
                    diffusedAccumulatedMPEPT=diffusedAccumulatedMPEPT+diffusedMPEPT;
                end
            end
            gridLB(ix,iy)=gridLB(ix,iy)-diffusedAccumulatedMPEPT;
        end
    end

    %Folic Acid Diffusion
    for ix = 1:j
        for iy=1:k
            diffusedAccumulatedFA=0;
            if(gridTempST(ix,iy)==0)
                continue;
            end
            if(isIndexValid(gridST,ix,iy-1))
                if(gridTempST(ix,iy-1)<gridTempST(ix,iy))
                    diffusedFA=diffusion1(gridTempST(ix,iy-1),gridTempST(ix,iy),0.01,0.1);
                    gridST(ix,iy-1)=gridST(ix,iy-1)+diffusedFA;
                    diffusedAccumulatedFA=diffusedAccumulatedFA+diffusedFA;
                end
            end

            if(isIndexValid(gridST,ix,iy+1))
                if(gridTempST(ix,iy+1)<gridTempST(ix,iy))
                    diffusedFA=diffusion1(gridTempST(ix,iy+1),gridTempST(ix,iy),0.01,0.1);
                    gridST(ix,iy+1)=gridST(ix,iy+1)+diffusedFA;
                    diffusedAccumulatedFA=diffusedAccumulatedFA+diffusedFA;
                end
            end

            if(isIndexValid(gridST,ix-1,iy))
                if(gridTempST(ix-1,iy)<gridTempST(ix,iy))
                    diffusedFA=diffusion1(gridTempST(ix-1,iy),gridTempST(ix,iy),0.01,0.1);
                    gridST(ix-1,iy)=gridST(ix-1,iy)+diffusedFA;
                    diffusedAccumulatedFA=diffusedAccumulatedFA+diffusedFA;
                end
            end

            if(isIndexValid(gridST,ix+1,iy))
                if(gridTempST(ix+1,iy)<gridTempST(ix,iy))
                    diffusedFA=diffusion1(gridTempST(ix+1,iy),gridTempST(ix,iy),0.01,0.1);
                    gridST(ix+1,iy)=gridST(ix+1,iy)+diffusedFA;
                    diffusedAccumulatedFA=diffusedAccumulatedFA+diffusedFA;
                end
            end
            gridST(ix,iy)=gridST(ix,iy)-diffusedAccumulatedFA;
        end
    end
folicAcidCF = gridST(idxR,idyR);
mpeptCF = gridLB(idxL,idyL);
end
figure('Name','ST Toy Model Concentrations','NumberTitle','off')
subplot(1,2,1);
plot(cell2mat(BiomassConcentrationsST))
title('Biomass')
xlabel('Timesteps') 
ylabel('Biomass concentrations (units)') 
subplot(1,2,2);
a1 = plot(cell2mat(metaboliteConcentrationsST(1,:)));b1="Lactose";
hold on
a2 = plot(cell2mat(metaboliteConcentrationsST(2,:)),'LineWidth',2,'LineStyle','--');b2="Folic Acid";
hold on
a3 = plot(cell2mat(metaboliteConcentrationsST(3,:)));b3="Other Products";
hold on
a4 = plot(cell2mat(metaboliteConcentrationsST(end,:)));b4="Peptides";
legend([a1;a2;a3;a4],b1,b2,b3,b4)
title('ST Concentrations');
xlabel('Timesteps') 
ylabel('Metabolite concentrations (units)') 
hold off

figure('Name','LB Toy Model Concentrations','NumberTitle','off')
subplot(1,2,1);
plot(cell2mat(BiomassConcentrationsLB))
title('Biomass')
xlabel('Timesteps') 
ylabel('Biomass concentrations (units)') 
subplot(1,2,2);
a1 = plot(cell2mat(metaboliteConcentrationsLB(1,:)),'LineWidth',2,'LineStyle','--');b1="Peptides";
hold on
a2 = plot(cell2mat(metaboliteConcentrationsLB(2,:)));b2="Folic Acid";
hold on
a3 = plot(cell2mat(metaboliteConcentrationsLB(3,:)));b3="Lactose";
hold on
a4 = plot(cell2mat(metaboliteConcentrationsLB(end,:)));b4="Other Products";
legend([a1;a2;a3;a4],b1,b2,b3,b4)
title('LB Concentrations');
xlabel('Timesteps') 
ylabel('Metabolite concentrations (units)') 
hold off

function oValue = isIndexValid(iTable,iIndexX,iIndexY)
try
    if(isa(iTable,'cell'))
        iTable{iIndexX,iIndexY}; %#ok<*VUNUS>
    else
        iTable(iIndexX,iIndexY);
    end
    oValue = true;
catch %#ok<CTCH>
    oValue = false;
end
end

