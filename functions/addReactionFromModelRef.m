function model = addReactionFromModelRef(model, ids, modelRef, includeRules, excludeRxnsByRelaxingDirection)

if nargin<4; includeRules = 0; end
if nargin<5; excludeRxnsByRelaxingDirection = 1; end;

if allMetsInCBMPYFormat(model.mets) && allMetsInCOBRAFormat(modelRef.mets)
    modelRef = transformModelToCBMPYFormat(modelRef);
end
if allMetsInCBMPYFormat(modelRef.mets) && allMetsInCOBRAFormat(model.mets)
    modelRef = transformModelToCOBRAFormat(modelRef);
end

if includeRules
   if ~isfield(modelRef,'grRules') && isfield(modelRef,'rules')
      modelRef = createGrRulesFromRules(modelRef); 
   end
   if ~isfield(model,'grRules') && isfield(model,'rules')
      model = createGrRulesFromRules(model); 
   end  
   if isfield(modelRef,'grRules') && ~isfield(modelRef,'rules')
      modelRef = createRulesFromgrRules(modelRef); 
   end
   if isfield(model,'grRules') && ~isfield(model,'rules')
      model = createRulesFromgrRules(model); 
   end  
end

if ~isempty(strfind(model.mets{1},'[')) && isempty(strfind(model.mets{1},'['))
    modelRef.mets = regexprep(modelRef.mets,'_((?!_).)+$','[$1]');
end
for i = 1:length(ids)
    pos = find(strcmp(modelRef.rxns, ids{i}));
    
    if ~isempty(pos) && length(pos)==1
        
        if ~isempty(find(strcmp(model.rxns, ids{i})));
            warning('A reaction with the same ID was found. It will not be added')
            continue;
        end
        
        eq = getRxn_cobraFormat(modelRef, pos);
        [rxnFormulaWasFound, posRxn] = reactionFormulaInModel(model, eq{1},excludeRxnsByRelaxingDirection);
        if rxnFormulaWasFound
            fprintf(['We couln''t add the reaction "' ids{i} '" with equation "' eq{1} '"\n' ...
                'from the reference model because reaction "' model.rxns{posRxn} '" in model has the same equation\n'])
            warning('A reaction with the reaction equation was found. It will not be added')
            continue;
        end
        newMets = setdiff(modelRef.mets(find(modelRef.S(:,pos))),model.mets);
        model = addReaction(model, ids{i}, 'reactionFormula', eq{1});
        posInModel = find(strcmp(model.rxns, ids{i}));
        if isfield(modelRef,'rxnNames') && isfield(model,'rxnNames') && ~isempty(modelRef.rxnNames{pos})
            model.rxnNames{posInModel} = modelRef.rxnNames{pos};
        end
        if isfield(modelRef,'subSystems') && isfield(model,'subSystems') && ~isempty(modelRef.subSystems{pos})
            model.subSystems{posInModel} = modelRef.subSystems{pos};
        end
        if isfield(modelRef,'rxnKEGGID') && isfield(model,'rxnKEGGID') && ~isempty(modelRef.rxnKEGGID{pos})
            model.rxnKEGGID{posInModel} = modelRef.rxnKEGGID{pos};
        end
        if isfield(modelRef,'rxnBiGGID') && isfield(model,'rxnBiGGID') && ~isempty(modelRef.rxnBiGGID{pos})
            model.rxnBiGGID{posInModel} = modelRef.rxnBiGGID{pos};
        end
        if isfield(modelRef,'rxnECNumbers') && isfield(model,'rxnECNumbers') && ~isempty(modelRef.rxnECNumbers{pos})
            model.rxnECNumbers{posInModel} = modelRef.rxnECNumbers{pos};
        end
        if isfield(modelRef,'rxnMetaNetXID') && isfield(model,'rxnMetaNetXID') && ~isempty(modelRef.rxnMetaNetXID{pos})
            model.rxnMetaNetXID{posInModel} = modelRef.rxnMetaNetXID{pos};
        end
        if isfield(modelRef,'rxnSBOTerms') && isfield(model,'rxnSBOTerms') && ~isempty(modelRef.rxnSBOTerms{pos})
            model.rxnSBOTerms{posInModel} = modelRef.rxnSBOTerms{pos};
        end
        if isfield(modelRef,'rules') && isfield(model,'rules') && ~isempty(modelRef.rules{pos}) && includeRules
            model.rules{posInModel} = modelRef.rules{pos};
        end
        if isfield(modelRef,'grRules') && isfield(model,'grRules') && ~isempty(modelRef.grRules{pos}) && includeRules
            model.grRules{posInModel} = modelRef.grRules{pos};
        end
        
        
        posMets = find(model.S(:,posInModel));
        posMetsInRef = cell2mat(arrayfun(@(x)find(strcmp(x,modelRef.mets)), model.mets(posMets), 'UniformOutput',false));
        for j = 1:length(posMets);
            posMetInModel_j = posMets(j);
            posMetInRef_j = posMetsInRef(j);
            if isfield(modelRef,'metNames') && isfield(model,'metNames') && ismember(model.mets{posMets(j)},newMets)
                model.metNames{posMetInModel_j} = modelRef.metNames{posMetInRef_j};
            end
            if isfield(modelRef,'metCharges') && isfield(model,'metCharges') && ~isempty(modelRef.metCharges(posMetInRef_j)) && isempty(model.metCharges(posMetInModel_j))
                model.metCharges{posMetInModel_j} = modelRef.metCharges{posMetInRef_j};
            end
            if isfield(modelRef,'metFormulas') && isfield(model,'metFormulas') && ~isempty(modelRef.metFormulas{posMetInRef_j}) && isempty(model.metFormulas{posMetInModel_j})
                model.metFormulas{posMetInModel_j} = modelRef.metFormulas{posMetInRef_j};
            end
            if isfield(modelRef,'metKEGGID') && isfield(model,'metKEGGID') && ~isempty(modelRef.metKEGGID{posMetInRef_j}) && isempty(model.metKEGGID{posMetInModel_j})
                model.metKEGGID{posMetInModel_j} = modelRef.metKEGGID{posMetInRef_j};
            end
            if isfield(modelRef,'metChEBIID') && isfield(model,'metChEBIID') && ~isempty(modelRef.metChEBIID{posMetInRef_j}) && isempty(model.metChEBIID{posMetInModel_j})
                model.metChEBIID{posMetInModel_j} = modelRef.metChEBIID{posMetInRef_j};
            end
            if isfield(modelRef,'metMetaNetXID') && isfield(model,'metMetaNetXID') && ~isempty(modelRef.metMetaNetXID{posMetInRef_j}) && isempty(model.metMetaNetXID{posMetInModel_j})
                model.metMetaNetXID{posMetInModel_j} = modelRef.metMetaNetXID{posMetInRef_j};
            end
            if isfield(modelRef,'metSBOTerms') && isfield(model,'metSBOTerms') && ~isempty(modelRef.metSBOTerms{posMetInRef_j}) && isempty(model.metSBOTerms{posMetInModel_j})
                model.metSBOTerms{posMetInModel_j} = modelRef.metSBOTerms{posMetInRef_j};
            end
        end
    else
        warning('ID was not found in modelRef database')
    end
end

if includeRules
   if isfield(model, 'grRules') && isfield(modelRef, 'grRules')
        model = createGenesFromGrRules(model);
   elseif isfield(model, 'rules') && isfield(modelRef, 'rules')
       model = createGrRulesFromRules(model); 
       model = createGenesFromGrRules(model);
   end
end

for i = 1:length(model.subSystems)
    if iscell(model.subSystems{i})
        model.subSystems{i} = model.subSystems{i}{1};
    end
    if isempty(model.subSystems{i})
        model.subSystems{i} = '';
    end
end

end