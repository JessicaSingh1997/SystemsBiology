function [growthRate,status] = maximizeGrowth(model,rxnID)
%function description
growthRate = "NA";
model = changeRxnBounds(model,rxnID,0.0001,'l');
fba = optimizeCbModel(model);
status = fba.stat;
if fba.f > 0.00000001
    growthRate = fba.f;
end
end
