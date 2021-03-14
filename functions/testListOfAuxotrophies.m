function dispResults = testListOfAuxotrophies
model = readSBML('/Users/jessica/SystemsBiology/S_thermophillus.xml',1000);
listOfReactions = model.rxns(find(findExcRxns(model)));
results = zeros(length(listOfReactions),1);
for i = 1:length(listOfReactions)
 result = testAuxotrophy(model, listOfReactions{i});
 results(i) = result;
end
dispResults = {listOfReactions num2cell(results)}
end