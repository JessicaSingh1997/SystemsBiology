function dispGrowthRate = testAuxotrophies
model = readSBML('/Users/jessica/SystemsBiology/S_thermophillus.xml',1000);
listOfReactions = model.rxns(find(findExcRxns(model)));
results = zeros(length(listOfReactions),1);
for i = 1:length(listOfReactions)
 result = testAuxotrophy(model, listOfReactions{i});
 results(i) = result;
 if result == 0
     dispGrowthRate = testAuxotrophicMetabolites(model,listOfReactions{i})
 end
end
end