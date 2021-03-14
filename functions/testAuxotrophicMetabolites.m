function modifiedGrowthRate = testAuxotrophicMetabolites(model,rxnID)
modifiedGrowthRate = 0;
model = changeRxnBounds(model,rxnID,-1001,'l');
fba = optimizeCbModel(model);
if fba.f > 0.6013
    modifiedGrowthRate = fba.f;
end
end