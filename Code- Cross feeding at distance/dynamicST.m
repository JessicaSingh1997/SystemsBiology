function[folicacidProduced,metaboliteConcentrationsST,BiomassConcentrationsST,mu] = dynamicST(n,n_rows, n_cols, n_metabolites_to_trackST,metaboliteConcentrationsST,BiomassConcentrationsST, model,dt,exchangeRxns_metabolites_to_trackST,caspepCF)

%Performs Dynamic FBA on ST and provides Folic acid for diffusion and recieves peptides from LB and adds to media
%
%[folicacidProduced,metaboliteConcentrationsST,BiomassConcentrationsST,mu] = dynamicFBAST(n,n_rows, n_cols, n_metabolites_to_trackST,metaboliteConcentrationsST,BiomassConcentrationsST, model,dt,exchangeRxns_metabolites_to_trackST,caspepCF)
%
%INPUT:
%   n - Number of iterations
%   n_metabolites_to_trackST - Number of exchange reactions being tracked
%   metaboliteConcentrationsST - initial metabolite concentrations 
%   BiomassConcentrationsST - initial biomass concentration
%   model - ST model used
%   dt - time step
%   exchangeRxns_metabolites_to_trackST - names of reactions being tracked
%   caspepCF - cross-fed casein peptide concentration
%
%OUTPUTS:
% folicacidProduced - folic acid produced for cross-feeding 
% metaboliteConcentrationsST - updated concentrations
% BiomassConcentrationsST - updated biomass 
% mu - growth rate
%

    for i = 1:n_rows
        for j = 1:n_cols
            concentrationsST = zeros(n_metabolites_to_trackST,1);
            for m = 1:n_metabolites_to_trackST
                c = metaboliteConcentrationsST{m,n};
                concentrationsST(m) = c(i,j);
            end
            biomassConcentrationST = BiomassConcentrationsST{n}(i,j);
            disp('')
            if all(concentrationsST==0)
                
            else
                [concentrationMatrixST, excRxnNamesST, timeVec,...
                    biomassVec,mu] = dynamicFBAm(model, exchangeRxns_metabolites_to_trackST,...
                    concentrationsST, biomassConcentrationST, dt, 1,exchangeRxns_metabolites_to_trackST);
                concentrationMatrixST = full(concentrationMatrixST);
                BiomassConcentrationsST{n+1}(i,j) = biomassVec(end);
                for m = 1:n_metabolites_to_trackST
                    if ismember(exchangeRxns_metabolites_to_trackST(m),excRxnNamesST)
                        pos =  find(strcmp(excRxnNamesST,exchangeRxns_metabolites_to_trackST(m)));  
                        metaboliteConcentrationsST{m,n+1}(i,j) = concentrationMatrixST(pos,end);
                    end 
                end
                pos_caspep =  find(strcmp(excRxnNamesST,exchangeRxns_metabolites_to_trackST(3)));
                concentrationMatrixST(pos_caspep,end) = concentrationMatrixST(pos_caspep,end)+ caspepCF; %cross-fed caspep added to concentrations
                metaboliteConcentrationsST{3,n+1}(i,j) = concentrationMatrixST(pos_caspep,end);
            end            
        end
        folicacidProduced=metaboliteConcentrationsST{7,n+1};%pos=18of37
    end
end