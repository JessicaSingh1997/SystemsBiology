function growthResults = testListOfGrowth
model = readSBML('/Users/jessica/SystemsBiology/LBUL.xml',1000);
listOfReactions = model.rxns(find(findExcRxns(model)));
results = zeros(length(listOfReactions),1);
statuses = zeros(length(listOfReactions),1);
for i = 1:length(listOfReactions)
 [result_GrowthRate,result_Status] = maximizeGrowth(model, listOfReactions{i});
 results(i) = result_GrowthRate;
 statuses(i) = result_Status;
end
growthResults = {listOfReactions num2cell(results) statuses};
end