function [diffusionFlux] = diffusion(concentrations1,concentrations2,D,X)

%D - Diffusion constant value - Can set to 0.0001 cm^2/s
%Idea - Link to table with D values for different metabolites
%D = 0.001;

%x - Distance between 2 cells of organism
%Idea - To define grid space - can be possible for mutiple cells to be 
%present in the same grid. Alternative - make cross section of grid area
%very small in order to reduce the number of dimensions. 
%X = 0.1;

%diffusionFlux - gives concentration flux of 
diffusionFlux = ((concentrations2 - concentrations1)* D)/X;