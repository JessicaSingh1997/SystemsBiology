function[model_LB] = LB
%Code for LB.mat file
%Changes are made to LBUL.xml file so that it replicates yogurt medium 

modelLB = readSBML('/Users/jessica/SystemsBiology/LBUL.xml',1000);

%To add CASPEP exchange reaction
modelLB = addReaction(modelLB, 'EX_caspep_e', 'metaboliteList', {'caspep[e]'},...
'stoichCoeffList', [-1]);

%Changes made to ensure Lac_D production
modelLB.lb(423)=0; %Makes DHORDi lb=0,lac_D = 0.2605
modelLB.ub(515)=0; %Makes POX2 ub=0, we get lac_D result but no succ,acetd
modelLB.lb(731)=0.001; %Force caspep production
modelLB.lb(720)=-0.002; %Cas uptake rate
%model.ub(721)=10;model.lb(721)=0.1; %CASPROTL 
%model.lb(731)=0; model.ub(731)=100; %EX_caspep_e

%Changes to make cas the main source of AA - Makes AA uptake 0. 323,340,375
modelLB.lb(362)=0;modelLB.lb(365)=0;modelLB.lb(370)=0;modelLB.lb(376)=0;%model.lb(375)=0;
modelLB.lb(377 )=0;modelLB.lb(379)=0; modelLB.lb(345)=0;modelLB.lb(347)=0;modelLB.lb(351)=0;
modelLB.lb(352)=0;modelLB.lb(355)=0;modelLB.lb(342)=0;modelLB.lb(341)=0;%model.lb(340)=0;
modelLB.lb(332)=0;modelLB.lb(324)=0;modelLB.lb(322)=0;modelLB.lb(320)=0;%model.lb(323)=0; 
model_LB=modelLB;
end




