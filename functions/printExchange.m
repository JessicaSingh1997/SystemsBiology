function printExchange(model)
excRxns = findExcRxns(model);
rxnText = printRxnFormula(model, model.rxns, false);
lb = model.lb(excRxns);
ub = model.ub(excRxns);
for k = 1:length(rxnText)
    fprintf('%s  :\tlb = %d\tub = %d\n', rxnText{k}, lb(k), ub(k));
end

