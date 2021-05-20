function [BiomassConcentrationsST,BiomassConcentrationsLB,metaboliteConcentrationsST,metaboliteConcentrationsLB,muST,muLB] = CrossFeeding(x,y,idyR,idyL,totalSimulationTime, dt)
%Performs cross-feeding between ST and LB
%
%USAGE: [BiomassConcentrationsST,BiomassConcentrationsLB,metaboliteConcentrationsST,metaboliteConcentrationsLB,muST,muLB] = CrossFeeding(x,y,idyR,idyL, totalSimulationTime, dt)
%
%INPUTS:
%       x                   - Number of rows in grid (height) in μm
%       y                   - Number of columns in grid (width) in μm
%       idyL                - Position of ST biomass on y axis (idxR
%                             (position of ST biomass on x axis - Has been
%                             set to middle row of grid))
%       idyR                - Position of LB biomass on y axis 
%                             (idxL - position of LB biomass on x axis - 
%                               Has been set to middle row of grid)
%       totalSimulationTime - Total time needed for simulation (Hours)
%       dt                  - Time step
%
%OUTPUTS:
%       BiomassConcentrationST      - Biomass of ST over time
%       BiomassConcentrationLB      - Biomass of LB over time
%       metaboliteConcentrationST   - Concentration of metabolites consumed/produced by of ST over time
%       metaboliteConcentrationLB   - Concentration of metabolites consumed/produced by of of LB over time
%       muST                        - Growth rate of ST
%       muLB                        - Growth rate of LB
%
%NOTE: Additional commands to plot graphs have been commented at end of code 

%Set up Grid Space - GridR/L-final concentrations GridTempR/L where changes are tracked

gridLB = zeros(x,y); %x and y coordinate of the surface where the LB is - Right
gridTempLB=zeros(x,y); %where changes in concentration are stored for LB (Right)
gridST = zeros(x,y); %x and y coordinate of the surface where the ST is - Left
gridTempST=zeros(x,y); %where changes in concentration are stored for ST (Left)

%Initialize the position of Biomass on Grid
%Setting distance between cells
[j,k]= size(gridLB);
%Postion of LB biomass at t=0
idxR=ceil(j/2);
%idyR=41;%ceil(3*k/4)(8-15);
%Position of ST biomass at t=0
idxL=ceil(j/2);
%idyL=26;%ceil(k/4)1-7;
gridDBR = gridLB;
gridDBL = gridST;
gridST(idxL,idyL)=0; %Initialize concentration of CF metabolite - Folic Acid to zero
gridLB(idxR,idyR)=0; %Initialize concentration of CF metabolite - mpept to zero

%Number of timesteps
%dt = 0.1; %hr
%totalSimulationTime = 30;
n_iterations = totalSimulationTime/dt;

%%%_______________________________ MODEL ST _______________________________%%%
modelST=ST;
%load('Model_LB') - load model from ST file, save as mat file to workspace to avoid repeated
%recomputation
modelST.lb(282) = -12; %Reduce lactose uptake
fbaST = optimizeCbModel(modelST);
mediaST = getMediaFromModel(modelST);
nutrientsST = regexprep(mediaST.reactions,{'EX_','_e'},{'',''});
positions_exchange_rxnsST = find(startsWith(modelST.rxns,'EX_'));
pos_exchange_with_fluxST = intersect(find(fbaST.x~=0), positions_exchange_rxnsST);
exchange_rxns_with_fluxST = modelST.rxns(pos_exchange_with_fluxST);
productsST =  regexprep(setdiff(exchange_rxns_with_fluxST, mediaST.reactions),{'EX_','_e'},{'',''});
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

initialbiomassConcentrationST= 0.005; %Initialize ST biomass
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

%%%________________________________________MODEL LB_________________________________%%%
modelLB=LB;
%load('Model_LB')%Generated from ST.m file with changes made to ST model
modelLB.lb(336)=-0.1;%Folate
fbaLB = optimizeCbModel(modelLB);
mediaLB = getMediaFromModel(modelLB);
nutrientsLB = regexprep(mediaLB.reactions,{'EX_','_e'},{'',''});
positions_exchange_rxnsLB = find(startsWith(modelLB.rxns,'EX_'));
pos_exchange_with_fluxLB = intersect(find(fbaLB.x~=0), positions_exchange_rxnsLB);
exchange_rxns_with_fluxLB = modelLB.rxns(pos_exchange_with_fluxLB);
productsLB =  regexprep(setdiff(exchange_rxns_with_fluxLB, mediaLB.reactions),{'EX_','_e'},{'',''});
metabolites_to_trackLB = union(nutrientsLB,productsLB);
n_metabolites_to_trackLB = length(metabolites_to_trackLB);
exchangeRxns_metabolites_to_trackLB = strcat('EX_', metabolites_to_trackLB, '_e');
grid0LB = zeros(n_rows,n_cols);

initialConcentrationsLB = cell(size(metabolites_to_trackLB));
for i = 1:length(metabolites_to_trackLB)
    initialConcentrationsLB{i} = zeros(n_rows,n_cols);
end

initialbiomassConcentrationLB= 0.008; %LB biomass gDW
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

%%%%%____________Nutrient Concentrations in ST Media________________%%%%%

for p = 1:length(metabolites_to_trackST)
    if ismember(metabolites_to_trackST(p),nutrientsST)
        initialConcentrationsST{3} = 0.0001; %CASPEP concentration
        initialConcentrationsST{19}= 100; %LCTS concentration
        initialConcentrationsST{4}=20; %Cit concentration
        initialConcentrationsST{15}=50; %Water concentration
        initialConcentrationsST{25}= 4; %Pi concentration
        initialConcentrationsST{p}= 1; %Other media components
        if p==34
            initialConcentrationsST{34}=10; %Urea concentration
        end      
    else
        initialConcentrationsST{p}= 0;
    end
end
for p = 1:n_metabolites_to_trackST
    metaboliteConcentrationsST{p,1} = initialConcentrationsST{p};
end

%%%%%____________Nutrient Concentrations in LB Media________________%%%%%

for q = 1:length(metabolites_to_trackLB)
    if ismember(metabolites_to_trackLB(q),nutrientsLB)
        initialConcentrationsLB{19}= 0; %Folic Acid Concentration
        initialConcentrationsLB{31}= 100; %LCTS concentration
        initialConcentrationsLB{40}= 20; %Glu
        initialConcentrationsLB{26}= 100; %Water
        initialConcentrationsLB{29}= 1; %
        initialConcentrationsLB{40}= 0.5; %NH4
        initialConcentrationsLB{42}= 8; %pi
        initialConcentrationsLB{q}= 0.1; %Other media components
        
    else
        initialConcentrationsLB{q} = 0;
    end
end
for q = 1:n_metabolites_to_trackLB
    metaboliteConcentrationsLB{q,1} = initialConcentrationsLB{q};
end


%%%%_______________________Cross Feeding_______________________%%%%

%Initialize metabolites being cross-fed
diffusedCASPEP=0;
diffusedAccumulatedCASPEP=0;
diffusedFA=0;
diffusedAccumulatedFA=0;
folicAcidCF=0;
caspepCF=0;
folicAcidCFlist=[];
muST=[]; %Returns growth rate of ST
muLB=[]; %Returns growth rate of LB


for i = 1:n_iterations
    diffusedAccumulatedCASPEP=0;
    diffusedAccumulatedFA=0;
    gridTempLB=gridLB;
    gridTempST=gridST;
    
    %Update LB Media concentration when diffused Folic Acid is available
    for q = 1:length(metabolites_to_trackLB)
        if ismember(metabolites_to_trackLB(q),nutrientsLB)
            initialConcentrationsLB{19}= folicAcidCF; %Cross-Fed Folic Acid
            initialConcentrationsLB{31}= 1; %LCTS concentration
            initialConcentrationsLB{q}= 0.1; %Other media components
        else
            initialConcentrationsLB{q} = 0;
        end
    end
    for q = 1:n_metabolites_to_trackLB
        if ismember(metabolites_to_trackLB(q),nutrientsLB)
            metaboliteConcentrationsLB{19,i} = initialConcentrationsLB{19};
        else
            metaboliteConcentrationsLB{q,1} = initialConcentrationsLB{q};
        end
    end
    
    
    %%_____Dynamic FBA of ST and LB__________%%
    
    if i~=n_iterations
        [folicacidProduced,metaboliteConcentrationsST,BiomassConcentrationsST,muST(end+1)] = dynamicST(i,n_rows, n_cols, n_metabolites_to_trackST,metaboliteConcentrationsST,BiomassConcentrationsST,modelST,dt,exchangeRxns_metabolites_to_trackST,caspepCF);
        [caspepProduced,metaboliteConcentrationsLB,BiomassConcentrationsLB,muLB(end+1)] = dynamicLB(i,n_rows, n_cols, n_metabolites_to_trackLB,metaboliteConcentrationsLB,BiomassConcentrationsLB, modelLB,dt,exchangeRxns_metabolites_to_trackLB);
    end
    
    gridST(idxL,idyL)=gridST(idxL,idyL)+folicacidProduced; %ST position of produced Folic Acid
    gridLB(idxR,idyR)=gridLB(idxR,idyR)+caspepProduced;%LB position of produced peptides (caspep)
    
    %Peptide diffusion in LB grid
    %Checks concentrations in all cells around it and distributes
    %concentrations evenly. 
    for ix = 1:j
        for iy=1:k
            countP=1;
            
            if isIndexValid(gridLB,ix,iy-1)
                countP=countP+1;
            end
            if isIndexValid(gridLB,ix,iy+1)
                countP=countP+1;
            end
            if isIndexValid(gridLB,ix-1,iy)
                countP=countP+1;
            end
            if isIndexValid(gridLB,ix+1,iy)
                countP=countP+1;
            end
            
            diffusedAccumulatedCASPEP=0;
            maxDiffusedP = gridTempLB(ix,iy)/countP;
            %Maximum amount of peptide that can be diffused in all 4
            %directions equally
            %Diffusion function called
            
            if(gridTempLB(ix,iy)==0)
                continue;
            end
            if(isIndexValid(gridLB,ix,iy-1))
                if(gridTempLB(ix,iy-1)<gridTempLB(ix,iy))
                    diffusedCASPEP=diffusion(gridTempLB(ix,iy-1),gridTempLB(ix,iy),5300000,1);
                    if diffusedCASPEP<maxDiffusedP
                        gridLB(ix,iy-1)=gridLB(ix,iy-1)+diffusedCASPEP;
                        diffusedAccumulatedCASPEP=diffusedAccumulatedCASPEP+diffusedCASPEP;
                    else
                        gridLB(ix,iy-1)=gridLB(ix,iy-1)+maxDiffusedP;
                        diffusedAccumulatedCASPEP=diffusedAccumulatedCASPEP+maxDiffusedP;
                    end
                end
            end
            
            if(isIndexValid(gridLB,ix,iy+1))
                if(gridTempLB(ix,iy+1)<gridTempLB(ix,iy))
                    diffusedCASPEP=diffusion(gridTempLB(ix,iy+1),gridTempLB(ix,iy),5300000,1);
                    if diffusedCASPEP<maxDiffusedP
                        gridLB(ix,iy+1)=gridLB(ix,iy+1)+diffusedCASPEP;
                        diffusedAccumulatedCASPEP=diffusedAccumulatedCASPEP+diffusedCASPEP;
                    else
                        gridLB(ix,iy+1)=gridLB(ix,iy+1)+maxDiffusedP;
                        diffusedAccumulatedCASPEP=diffusedAccumulatedCASPEP+maxDiffusedP;
                    end
                end
            end
            
            if(isIndexValid(gridLB,ix-1,iy))
                if(gridTempLB(ix-1,iy)<gridTempLB(ix,iy))
                    diffusedCASPEP=diffusion(gridTempLB(ix-1,iy),gridTempLB(ix,iy),5300000,1);
                    if diffusedCASPEP<maxDiffusedP
                        gridLB(ix-1,iy)=gridLB(ix-1,iy)+diffusedCASPEP;
                        diffusedAccumulatedCASPEP=diffusedAccumulatedCASPEP+diffusedCASPEP;
                    else
                        gridLB(ix-1,iy)=gridLB(ix-1,iy)+maxDiffusedP;
                        diffusedAccumulatedCASPEP=diffusedAccumulatedCASPEP+maxDiffusedP;
                    end
                end
            end
            
            if(isIndexValid(gridLB,ix+1,iy))
                if(gridTempLB(ix+1,iy)<gridTempLB(ix,iy))
                    diffusedCASPEP=diffusion(gridTempLB(ix+1,iy),gridTempLB(ix,iy),5300000,1);
                    if diffusedCASPEP<maxDiffusedP
                        gridLB(ix+1,iy)=gridLB(ix+1,iy)+diffusedCASPEP;
                        diffusedAccumulatedCASPEP=diffusedAccumulatedCASPEP+diffusedCASPEP;
                    else
                        gridLB(ix+1,iy)=gridLB(ix+1,iy)+maxDiffusedP;
                        diffusedAccumulatedCASPEP=diffusedAccumulatedCASPEP+maxDiffusedP;
                    end
                end
            end
            gridLB(ix,iy)=gridLB(ix,iy)-diffusedAccumulatedCASPEP; %Update current cell in LB gird
        end
    end
    
    %Folic Acid Diffusion
    for ix = 1:j
        for iy=1:k
            countF=1;
            
            if isIndexValid(gridST,ix,iy-1)
                countF=countF+1;
            end
            if isIndexValid(gridST,ix,iy+1)
                countF=countF+1;
            end
            if isIndexValid(gridST,ix-1,iy)
                countF=countF+1;
            end
            if isIndexValid(gridST,ix+1,iy)
                countF=countF+1;
            end
            diffusedAccumulatedFA=0;
            maxDiffusedF = gridTempST(ix,iy)/countF;
            
            if(gridTempST(ix,iy)==0)
                continue;
            end
            if(isIndexValid(gridST,ix,iy-1))
                if(gridTempST(ix,iy-1)<gridTempST(ix,iy))
                    diffusedFA=diffusion(gridTempST(ix,iy-1),gridTempST(ix,iy),1800000,1);
                    if diffusedFA<maxDiffusedF
                        gridST(ix,iy-1)=gridST(ix,iy-1)+diffusedFA;
                        diffusedAccumulatedFA=diffusedAccumulatedFA+diffusedFA;
                    else
                        gridST(ix,iy-1)=gridST(ix,iy-1)+maxDiffusedF;
                        diffusedAccumulatedFA=diffusedAccumulatedFA+maxDiffusedF;
                    end
                end
            end
            
            if(isIndexValid(gridST,ix,iy+1))
                if(gridTempST(ix,iy+1)<gridTempST(ix,iy))
                    diffusedFA=diffusion(gridTempST(ix,iy+1),gridTempST(ix,iy),1800000,1);
                    if diffusedFA<maxDiffusedF
                        gridST(ix,iy+1)=gridST(ix,iy+1)+diffusedFA;
                        diffusedAccumulatedFA=diffusedAccumulatedFA+diffusedFA;
                    else
                        gridST(ix,iy+1)=gridST(ix,iy+1)+maxDiffusedF;
                        diffusedAccumulatedFA=diffusedAccumulatedFA+maxDiffusedF;
                    end
                end
            end
            
            if(isIndexValid(gridST,ix-1,iy))
                if(gridTempST(ix-1,iy)<gridTempST(ix,iy))
                    diffusedFA=diffusion(gridTempST(ix-1,iy),gridTempST(ix,iy),1800000,1);
                    if diffusedFA<maxDiffusedF
                        gridST(ix-1,iy)=gridST(ix-1,iy)+diffusedFA;
                        diffusedAccumulatedFA=diffusedAccumulatedFA+diffusedFA;
                    else
                        gridST(ix-1,iy)=gridST(ix-1,iy)+maxDiffusedF;
                        diffusedAccumulatedFA=diffusedAccumulatedFA+maxDiffusedF;
                    end
                end
            end
            
            if(isIndexValid(gridST,ix+1,iy))
                if(gridTempST(ix+1,iy)<gridTempST(ix,iy))
                    diffusedFA=diffusion(gridTempST(ix+1,iy),gridTempST(ix,iy),1800000,1);
                    if diffusedFA<maxDiffusedF
                        gridST(ix+1,iy)=gridST(ix+1,iy)+diffusedFA;
                        diffusedAccumulatedFA=diffusedAccumulatedFA+diffusedFA;
                    else
                        gridST(ix+1,iy)=gridST(ix+1,iy)+maxDiffusedF;
                        diffusedAccumulatedFA=diffusedAccumulatedFA+maxDiffusedF;
                    end
                end
            end
            gridST(ix,iy)=gridST(ix,iy)-diffusedAccumulatedFA; %Update current cell in ST list
        end
    end
    
    folicAcidCF = gridST(idxR,idyR); %Folic acid cross fed to LB, concentration of Folic acid on grid element where LB biomass is defined
    caspepCF = gridLB(idxL,idyL); %Peptides cross fed to ST, concentration of Folic acid on grid element where ST biomass is defined
    
end


%{
%%%% CODE FOR PLOTS
1. ST Tracked metabolite concentrations
figure('Name','S. thermophillus Concentrations','NumberTitle','off')
subplot(1,2,1);
plot(cell2mat(BiomassConcentrationsST))
title('Biomass')
xlabel('Timesteps')
ylabel('Biomass concentrations of ST(units)')
subplot(1,2,2);
a1 = plot(cell2mat(metaboliteConcentrationsST(8,:)));b1="Formate";
hold on
a2 = plot(cell2mat(metaboliteConcentrationsST(7,:)),'LineWidth',2,'LineStyle','--');b2="Folic Acid";
hold on
a3 = plot(cell2mat(metaboliteConcentrationsST(5,:)));b3="Carbon Dioxide";
hold on
a4 = plot(cell2mat(metaboliteConcentrationsST(3,:)));b4="Casein Peptides";
hold on
a5 = plot(cell2mat(metaboliteConcentrationsST(19,:)));b5="Lactose";
hold on
a6 = plot(cell2mat(metaboliteConcentrationsST(18,:)));b6="Lactic Acid";
hold on
a7 = plot(cell2mat(metaboliteConcentrationsST(4,:)));b7="Cit";
hold on
a8 = plot(cell2mat(metaboliteConcentrationsST(15,:)),'LineWidth',2,'LineStyle','--');b8="H20";
hold on
a9 = plot(cell2mat(metaboliteConcentrationsST(25,:)));b9="pi";
legend([a1;a2;a3;a4;a5;a6;a7;a8;a9],b1,b2,b3,b4,b5,b6,b7,b8,b9)
title('ST Concentrations');
xlabel('Timesteps')
ylabel('Metabolite concentrations (units)')
hold off


2. LB Tracked Metabolite concentrations 
figure('Name','L. bulgaricus Concentrations','NumberTitle','off')
subplot(1,2,1);
plot(cell2mat(BiomassConcentrationsLB))
title('Biomass')
xlabel('Timesteps')
ylabel('Biomass concentrations of LB(units)')
subplot(1,2,2);
a1 = plot(cell2mat(metaboliteConcentrationsLB(10,:)),'LineWidth',2,'LineStyle','--');b1="Casein Peptides";
hold on
a2 = plot(cell2mat(metaboliteConcentrationsLB(19,:)));b2="Folic Acid";
hold on
a3 = plot(cell2mat(metaboliteConcentrationsLB(30,:)));b3="Lactic acid";
hold on
a4 = plot(cell2mat(metaboliteConcentrationsLB(31,:)));b4="Lactose";
legend([a1;a2;a3;a4],b1,b2,b3,b4)
title('LB Concentrations');
xlabel('Timesteps')
ylabel('Metabolite concentrations (units)')
hold off


3. ST and LB Biomass
figure('Name','Biomass Concentrations','NumberTitle','off')
a1 = plot(cell2mat(BiomassConcentrationsST));b1="ST";
yyaxis right
hold on
a2 = plot(cell2mat(BiomassConcentrationsLB));b2="LB";
legend([a1;a2],b1,b2)
hold off


4. Heatmaps for diffusion
%heatmap(gridTempLB, "ColorScaling", "log",'Colormap',autumn),h.Title = 'Diffusion of Casein Peptides from LB';
%heatmap(gridTempST, "ColorScaling", "log",'Colormap',winter),h.Title = 'Diffusion of Folic Acid from ST';


5. Plot Amino Acids
figure('Name','Amino Acid concentrations','NumberTitle','off')
a1 = plot(cell2mat(metaboliteConcentrationsST(1,:)));b1="Ala";
hold on
a2 = plot(cell2mat(metaboliteConcentrationsST(2,:)),'LineWidth',2,'LineStyle','--');b2="Arg";
hold on
a3 = plot(cell2mat(metaboliteConcentrationsST(6,:)));b3="Cys";
hold on
a4 = plot(cell2mat(metaboliteConcentrationsST(11,:)));b4="Glu";
hold on
a5 = plot(cell2mat(metaboliteConcentrationsST(12,:)));b5="Gly";
hold on
a6 = plot(cell2mat(metaboliteConcentrationsST(16,:)));b6="His";
hold on
a7 = plot(cell2mat(metaboliteConcentrationsST(17,:)));b7="Ile";
hold on
a8 = plot(cell2mat(metaboliteConcentrationsST(20,:)),'LineWidth',2,'LineStyle','--');b8="Leu";
hold on
a9 = plot(cell2mat(metaboliteConcentrationsST(21,:)));b9="Lys";
hold on
a10 = plot(cell2mat(metaboliteConcentrationsST(22,:)));b10="Met";
hold on
a11 = plot(cell2mat(metaboliteConcentrationsST(24,:)));b11="Phe";
hold on
a12 = plot(cell2mat(metaboliteConcentrationsST(27,:)));b12="Pro";
hold on
a13 = plot(cell2mat(metaboliteConcentrationsST(32,:)),'LineWidth',2,'LineStyle','--');b13="Trp";
hold on
a14 = plot(cell2mat(metaboliteConcentrationsST(33,:)),'LineWidth',2,'LineStyle','--');b14="Tyr";
hold on
a15 = plot(cell2mat(metaboliteConcentrationsST(35,:)),'LineWidth',2,'LineStyle','--');b15="Val";
legend([a1;a2;a3;a4;a5;a6;a7;a8;a9;a10;a11;a12;a13;a14;a15],b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b13,b14,b15)

%}

BiomassConcentrationsST = cell2mat(BiomassConcentrationsST);
BiomassConcentrationsLB = cell2mat(BiomassConcentrationsLB);
metaboliteConcentrationsST=cell2mat(metaboliteConcentrationsST);
metaboliteConcentrationsLB=cell2mat(metaboliteConcentrationsLB);

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
end

