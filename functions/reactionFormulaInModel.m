function [rxnFormulaWasFound, posRxn, rxnOppositeDirection, matchWithoutRemovingProtons, ...
    additionalMatchsIfRemovingProtons, posRxnsFullMatch, posRxnsPartialMatch, validPartialMatchs] ...
    = reactionFormulaInModel(model, equation, allowProductsByReactants, allowProtonMismatch)

if nargin < 4
    allowProtonMismatch = 0;
end

rxnFormulaWasFound = 0;
posRxn = [];
posRxnsFullMatch = [];
posRxnsPartialMatch = [];
validPartialMatchs = [];
rxnOppositeDirection = 0;
matchWithoutRemovingProtons = 0;
additionalMatchsIfRemovingProtons = 0;

if ischar(equation)
    [metaboliteList, stoichCoeffList, ~] = parseRxnFormula(equation);
elseif iscell(equation)
    [metaboliteList, stoichCoeffList, ~] = parseRxnFormula(equation{1});
    equation = equation{1};
elseif isnumeric(equation)
    posMets = find(model.S(:,equation));
    metaboliteList = model.mets(posMets);
    stoichCoeffList = full(model.S(posMets, equation));
    equation = getRxn_cobraFormat(model, equation);
    equation = equation{1};
elseif isstruct(equation)
    modelRef = equation.model;
    posRxnRef = find(strcmp(modelRef.rxns,equation.rxn));
    posMets = find(modelRef.S(:,posRxnRef));
    metaboliteList = modelRef.mets(posMets);
    stoichCoeffList = full(modelRef.S(posMets, posRxnRef));
end

[isInModel, posMets] = ismember(metaboliteList, model.mets);

if any(~isInModel)
    return
else
    
    if allowProtonMismatch && length(find(~cellfun(@isempty, regexp(metaboliteList, '^cpd00067_.*$|h_.*$|^PROTON_.*$|^C00080_.*$'))))~=length(posMets)
        posMetsWithProtons = posMets;
        stoichCoeffListWithProtons = stoichCoeffList;
        modelWithProtons = model;
        
        posProtons = find(~cellfun(@isempty, regexp(metaboliteList, '^cpd00067_.*$|h_.*$|^PROTON_.*$|^C00080_.*$')));
        kept = true(size(metaboliteList)); kept(posProtons) = false;
        posMets = posMets(kept);
        stoichCoeffList = stoichCoeffList(kept);
        metaboliteList = model.mets(posMets);
        posProtonsInModel = find(~cellfun(@isempty, regexp(model.mets, '^cpd00067_.*$|h_.*$|^PROTON_.*$|^C00080_.*$')));
        model = removeMetabolites(model, model.mets(posProtonsInModel),0);
        posMets = cell2mat(arrayfun(@(x)find(strcmp(x,model.mets)),metaboliteList,'UniformOutput',false))';
    end
    int = find(model.S(posMets(1),:));
    intersectedMets = 1;
    n_restantes = length(posMets) - intersectedMets;
    
    while n_restantes > 0 && ~isempty(int)
        intersectedMets = intersectedMets + 1;
        n_restantes = length(posMets) - intersectedMets;
        int = intersect(int, find(model.S(posMets(intersectedMets),:)));
    end
    if isempty(int)
        return;
    end
end
int = int(arrayfun(@(x) length(find(model.S(:,x))), int) == length(posMets));
posRxn = zeros(length(int),1);
rxnOppositeDirection = zeros(length(int),1);
cont = 0;

for i = 1:length(int)
    stoi = full(model.S(posMets,int(i)));
    if size(stoichCoeffList,1) ~= size(stoi,1); stoichCoeffList = stoichCoeffList'; end;
    if all(stoi == stoichCoeffList)
        cont = cont+1;
        posRxn(cont) = int(i);
    elseif (all(stoi == -stoichCoeffList) && allowProductsByReactants)
        cont = cont+1;
        rxnOppositeDirection(cont) = 1;
        posRxn(cont) = int(i);
    end
end

if cont > 0
    rxnFormulaWasFound = 1;
    posRxn = posRxn(1:cont);
    rxnOppositeDirection = rxnOppositeDirection(1:cont);
    
    if allowProtonMismatch && length(find(~cellfun(@isempty, regexp(metaboliteList, '^cpd00067_.*$|h_.*$|^PROTON_.*$|^C00080_.*$'))))~=length(posMets)
        
        int = find(modelWithProtons.S(posMetsWithProtons(1),:));
        intersectedMets = 1;
        n_restantes = length(posMetsWithProtons) - intersectedMets;
        
        while n_restantes > 0 && ~isempty(int)
            intersectedMets = intersectedMets + 1;
            n_restantes = length(posMetsWithProtons) - intersectedMets;
            int = intersect(int, find(modelWithProtons.S(posMetsWithProtons(intersectedMets),:)));
        end
        if ~isempty(int)
            int = int(arrayfun(@(x) length(find(modelWithProtons.S(:,x))), int) == length(posMetsWithProtons));
            posRxn2 = zeros(length(int),1);
            rxnOppositeDirection2 = zeros(length(int),1);
            cont2 = 0;
            
            for i = 1:length(int)
                stoi = full(modelWithProtons.S(posMetsWithProtons,int(i)));
                if size(stoichCoeffListWithProtons,1) ~= size(stoi,1); stoichCoeffListWithProtons = stoichCoeffListWithProtons'; end;
                if all(stoi == stoichCoeffListWithProtons)
                    cont2 = cont2+1;
                    posRxn2(cont2) = int(i);
                elseif (all(stoi == -stoichCoeffListWithProtons) && allowProductsByReactants)
                    cont2 = cont2+1;
                    rxnOppositeDirection2(cont2) = 1;
                    posRxn2(cont2) = int(i);
                end
            end
            if cont2>0
                posRxn2 = posRxn2(1:cont2);
                matchWithoutRemovingProtons = 1;
                posRxnsFullMatch = posRxn2;
                if ~isempty(setdiff(posRxn, posRxn2))
                    additionalMatchsIfRemovingProtons = 1;
                end
            end
            if cont>cont2
                if cont2==0
                    posRxnsPartialMatch = posRxn;
                else
                    posRxnsPartialMatch = setdiff(posRxn, posRxn2);
                end
            end
        else
            posRxnsPartialMatch = posRxn;
        end
        validPartialMatchs = cellfun(@(x) validateRxnsDifferentProtonStoichiometry(equation, x), getRxn_cobraFormat(modelWithProtons, posRxn));
    else
        posRxnsFullMatch = posRxn;
        matchWithoutRemovingProtons = 1;
        eqs = getRxn_cobraFormat(model, posRxn);
        validPartialMatchs = cellfun(@(x) validateRxnsDifferentProtonStoichiometry(equation, x), eqs);
    end
end 

end
