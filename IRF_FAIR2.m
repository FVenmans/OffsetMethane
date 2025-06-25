%This function creates a vector with the temperature impact response function of a project or emission impulse response. 
%INPUTS
%E_base: Emissions in the baseline
%E_project: Emissions with the project (E_base+1 in first period for standard IRF)
%F_ext : other gas forcing (for methane this could be approximated by CO2 forcing or total forcing)
%R0 5x1 vector with initial carbon boxes and cumulative emissions
%S0 3x1 initial temperature in temperature boxes (sum is initial
%temperature)
%paramR and paramT see FAIR2 function
% gamma: total damages are gamma/2*T²
% discount: discount rate for SCC
% g: growth rate of the economy
%OUTPUTS
% SVO: Social value of the project (negative for an absorption, positive for emission) as a proportion of GDP

function [T_IRF, F_IRF, Damage_IRF,PV_Damage_IRF, GWP, SVO]=IRF_FAIR2(E_base, E_project, F_ext, R0, S0, paramR, paramT,gamma,discount,g) 
t_end=size(E_base,2);
%without project
[T_base,F_base,~,~]=FAIR2(E_base, F_ext, R0, S0, paramR, paramT);
%with project
[T_project,F_project,~,~]=FAIR2(E_project, F_ext, R0, S0, paramR, paramT);

T_IRF=(T_project-T_base);
F_IRF=(F_project-F_base);
Damage_IRF=T_IRF.*gamma.*T_base;
PV_Damage_IRF=Damage_IRF.*exp(-(discount-g)*(0:t_end-1));
GWP=sum(F_IRF); 
SVO=sum(PV_Damage_IRF);

