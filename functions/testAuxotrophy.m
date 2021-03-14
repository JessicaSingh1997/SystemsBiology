function isAbleToGrow = testAuxotrophy(model,rxnID)
%function description
isAbleToGrow = 0;
model = changeRxnBounds(model,rxnID,0,'l');
fba = optimizeCbModel(model);
if fba.f > 0.001
    isAbleToGrow = 1;
end
end