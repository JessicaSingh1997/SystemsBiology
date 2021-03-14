function[mpeptProduced,metaboliteConcentrationsLB,BiomassConcentrationsLB] = dynamicLB(n,n_rows, n_cols, n_metabolites_to_trackLB,metaboliteConcentrationsLB,BiomassConcentrationsLB, modelLB,dt,exchangeRxns_metabolites_to_trackLB)
for i = 1:n_rows
        for j = 1:n_cols
            concentrationsLB = zeros(n_metabolites_to_trackLB,1);
            for m = 1:n_metabolites_to_trackLB
                c = metaboliteConcentrationsLB{m,n};
                concentrationsLB(m) = c(i,j);
            end
            biomassConcentrationLB = BiomassConcentrationsLB{n}(i,j);
            disp('')
            if all(concentrationsLB==0)
                
            else
                [concentrationMatrixLB, excRxnNamesLB, timeVec,...
                    biomassVec] = dynamicFBA(modelLB, exchangeRxns_metabolites_to_trackLB,...
                    concentrationsLB, biomassConcentrationLB, dt, 1, exchangeRxns_metabolites_to_trackLB);
                concentrationMatrixLB = full(concentrationMatrixLB);
                BiomassConcentrationsLB{n+1}(i,j) = biomassVec(end);
                for m = 1:n_metabolites_to_trackLB
                    if ismember(exchangeRxns_metabolites_to_trackLB(m),excRxnNamesLB)
                        pos =  find(strcmp(excRxnNamesLB,exchangeRxns_metabolites_to_trackLB(m)));
                        metaboliteConcentrationsLB{m,n+1}(i,j) = concentrationMatrixLB(pos,end);
                        %if folicAcidCF~=0
                            %concentrationMatrixLB(2,end) = concentrationMatrixLB(2,end)+ folicAcidCF;
                        %end
                    end
                end
            end               
        end
        mpeptProduced=metaboliteConcentrationsLB{1,n+1};
end
end

   
