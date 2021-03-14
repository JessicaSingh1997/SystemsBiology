function[folicacidProduced,metaboliteConcentrationsST,BiomassConcentrationsST] = dynamicST(n,n_rows, n_cols, n_metabolites_to_trackST,metaboliteConcentrationsST,BiomassConcentrationsST, model,dt,exchangeRxns_metabolites_to_trackST,mpeptCF)
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
                    biomassVec] = dynamicFBA(model, exchangeRxns_metabolites_to_trackST,...
                    concentrationsST, biomassConcentrationST, dt, 1,exchangeRxns_metabolites_to_trackST);
                concentrationMatrixST = full(concentrationMatrixST);
                BiomassConcentrationsST{n+1}(i,j) = biomassVec(end);
                for m = 1:n_metabolites_to_trackST
                    if ismember(exchangeRxns_metabolites_to_trackST(m),excRxnNamesST)
                        pos =  find(strcmp(excRxnNamesST,exchangeRxns_metabolites_to_trackST(m)));
                        metaboliteConcentrationsST{m,n+1}(i,j) = concentrationMatrixST(pos,end);
                        if mpeptCF~=0
                            concentrationMatrixST(4,end) = concentrationMatrixST(4,end)+ mpeptCF;
                        end
                    end
                end
            end            
        end
        folicacidProduced=metaboliteConcentrationsST{2,n+1};
    end
end

    