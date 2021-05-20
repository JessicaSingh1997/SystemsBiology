function[model_ST] = ST
%Code for ST.mat file
%Changes are made to S_thermophillus_v1_03.xml file so that it replicates
%yogurt medium 

modelST = readSBML('/Users/jessica/SystemsBiology/S_thermophillus_v1_03.xml',1000);
modelLB = readSBML('/Users/jessica/SystemsBiology/LBUL.xml',1000);

modelLB = addReaction(modelLB, 'EX_caspep_e', 'metaboliteList', {'caspep[e]'},...
'stoichCoeffList', [-1]);

%Remove MPEPT and replace with CASPEP as media
modelST = removeRxns(modelST,{'PEPTabc','CASL','CASabc','PEPTL','EX_mpept_e_','EX_caspep_e_'});
modelST = addReactionFromModelRef(modelST,{'CASL'},modelLB);
modelST = addReactionFromModelRef(modelST,{'CASabc'},modelLB);
modelST = addReactionFromModelRef(modelST,{'EX_caspep_e'},modelLB);

%Force production of Folic Acid and Lac_L
modelST.lb(518)=1;%Force production of lac_L
modelST.lb(14) =-1000; modelST.ub(14)= 1000; modelST.lb(312)=0.1; modelST.ub(312)=1000; %to force ST to produce fol - FOLt, EX_fol
model_ST=modelST;
end
