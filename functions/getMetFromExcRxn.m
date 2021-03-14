function [posMet, met] = getMetFromExcRxn(model, posExc)

if ischar(posExc)
    posExc = getPosOfElementsInArray({posExc}, model.rxns);    
end
posMet = find(model.S(:,posExc));
met = model.mets(posMet);

end