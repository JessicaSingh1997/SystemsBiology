function [posExcRxn, excRxns, posMets, excMets] = findExcRxnsWithIDs(model,option)

if nargin<2 || isempty(option); option = 1; end
if option ==1 % all
    posExcRxn = find(~cellfun(@isempty, regexp(model.rxns,'^EX_|^Ex_|^ex_|^R_EX_|^R_Ex_|^R_ex_')));
elseif option ==2 % just backward
    posExcRxn = find(~cellfun(@isempty, regexp(model.rxns,'^EX_.*_b$|^Ex_.*_b$|^ex_.*_b$|^R_EX_.*_b$|^R_Ex_.*_b$|^R_ex_.*_b$')));
elseif option ==3% just forward
    posExcRxn = find(~cellfun(@isempty, regexp(model.rxns,'^EX_.*_f$|^Ex_.*_f$|^ex_.*_f$|^R_EX_.*_f$|^R_Ex_.*_f$|^R_ex_.*_f$')));
elseif option ==4 %just backward, community modeling, species
    posExcRxn = find(~cellfun(@isempty, regexp(model.rxns,'^EX_.*_b_species.*$|^Ex_.*_b_species.*$|^ex_.*_b_species.*$|^R_EX_.*_b_species.*$|^R_Ex_.*_b_species.*$|^R_ex_.*_b_species.*$')));
elseif option ==5 %just forward, community modeling, species
    posExcRxn = find(~cellfun(@isempty, regexp(model.rxns,'^EX_.*_f_species.*$|^Ex_.*_f_species.*$|^ex_.*_f_species.*$|^R_EX_.*_f_species.*$|^R_Ex_.*_f_species.*$|^R_ex_.*_f_species.*$')));
elseif option ==6 %just uptake, community modeling, shared
    posExcRxn = find(~cellfun(@isempty, regexp(model.rxns,'^EX_i_shared_.*_e$|^Ex_i_shared_.*_e$|^ex_i_shared_.*_e$|^R_EX_i_shared_.*_e$|^R_Ex_i_shared_.*_e$|^R_ex_i_shared_.*_e$')));
elseif option ==7 %just secretion, community modeling, shared
    posExcRxn = find(~cellfun(@isempty, regexp(model.rxns,'^EX_o_shared_.*_e$|^Ex_o_shared_.*_e$|^ex_o_shared_.*_e$|^R_EX_o_shared_.*_e$|^R_Ex_o_shared_.*_e$|^R_ex_o_shared_.*_e$')));

end
excRxns = model.rxns(posExcRxn);
posMets = arrayfun(@(x) getMetFromExcRxn(model, x), posExcRxn);
excMets = model.mets(posMets);
end