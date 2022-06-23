function [p,flowrate,flowresult]=ferncodes_solverpressureNLFVPP(nflag,...
                                  parameter,kmap,wells,viscosity,V,Sw,N,...
                                  p_old,contnorm,wight,s)
%Define global parameters
global acel;

% interpolação nos nós ou faces
[pinterp]=ferncodes_pressureinterpNLFVPP(p_old,nflag,wight,s);
[M,I]=ferncodes_assemblematrixNLFVPP(pinterp,parameter,viscosity,contnorm);
%--------------------------------------------------------------------------
%Add a source therm to independent vector "mvector" 

%Often it may change the global matrix "M"
[M_old,RHS_old] = addsource(sparse(M),I,wells);
%% full Picard iteration
if strcmp(acel,'FPI')
    [p,flowrate,flowresult]=ferncodes_iterpicard(M_old,RHS_old,...
        parameter,wight,s,p_old,nflag,wells,viscosity);
elseif strcmp(acel,'AA')
    %% Picard-Anderson Acceleration
    [p,flowrate,flowresult]=ferncodes_iterpicardANLFVPP2(M_old,RHS_old,...
        parameter,wight,s,p_old,nflag,wells,viscosity,0,contnorm);
end
end