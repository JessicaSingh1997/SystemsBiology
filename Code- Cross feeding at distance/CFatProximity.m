%%% TO RUN MUTUAL INTERACTIONS ON ST AND LB BY CROSS-FEEDING %%%

%Select the position at whch ST and LB will be placed in grid to change distance between them - each element is 1μm wide. 
%Can enter mutiple positions to generate different distances
%idx - Position of ST
%idy - Positon of LB

idx=[24,25,26,27,28,29,30,44];
idy=[24,23,22,21,20,19,18,4];

%Distance is given by (idx(i) - idy(i))
%Distances in above case are 0,2,4,6,8,10,12,40 μm

BiomassConcentrationsST=[];
BiomassConcentrationsLB=[];
metaboliteConcentrationsST=[];
metaboliteConcentrationsLB=[];
muST=[];
muLB=[];

%From - Cross Feeding function 
%Performs cross-feeding between ST and LB
%
%USAGE: [BiomassConcentrationsST,BiomassConcentrationsLB,metaboliteConcentrationsST,metaboliteConcentrationsLB,muST,muLB] = CrossFeeding(x,y,idyR,idyL, totalSimulationTime, dt)
%
%INPUTS:
%       x                   - Number of rows in grid (height) in μm
%       y                   - Number of columns in grid (width) in μm
%       idyL                - Position of ST biomass on y axis (idxR
%                             (position of ST biomass on x axis - Has been
%                             set to middle row of grid))
%       idyR                - Position of LB biomass on y axis 
%                             (idxL - position of LB biomass on x axis - 
%                               Has been set to middle row of grid)
%       totalSimulationTime - Total time needed for simulation (Hours)
%       dt                  - Time step
%
%OUTPUTS:
%       BiomassConcentrationST      - Biomass of ST over time for multiple distances 
%       BiomassConcentrationLB      - Biomass of LB over time for multiple distances 
%       metaboliteConcentrationST   - Concentration of metabolites consumed/produced by of ST over time
%       metaboliteConcentrationLB   - Concentration of metabolites consumed/produced by of of LB over time
%       muST                        - Growth rate of ST for multiple distances 
%       muLB                        - Growth rate of LB for multiple distances 

for i = 1:size(idx,2)
    [BiomassConcentrationsST(end+1,:), BiomassConcentrationsLB(end+1,:),metaboliteConcentrationsST(end+1,:,:),metaboliteConcentrationsLB(end+1,:,:),muST(end+1,:),muLB(end+1,:)]= CrossFeeding(30,50,idx(i),idy(i),1,0.1);
    
end