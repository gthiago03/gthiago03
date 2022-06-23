
function [p,flowrate,flowresult]=ferncodes_solverpressureDMP(nflagface,...
    parameter,wells,mobility,weightDMP,p_old,w,s,nflagno,gravelem,gravpoint,gravface)
global acel gravresult gravrate keygravity P_old P_new A_old maxiter
% interpolação nos nós ou faces
[pinterp]=ferncodes_pressureinterpHP(p_old,nflagface,parameter,weightDMP);
% Montagem da matriz global

[M,I]=ferncodes_assemblematrixDMP(p_old,pinterp,0,parameter,weightDMP,mobility,gravelem,gravpoint,gravface);

%--------------------------------------------------------------------------
%Add a source therm to independent vector "mvector"

%Often it may change the global matrix "M"
[M_old,RHS_old] = addsource(sparse(M),I,wells);

%P_old = RHS_old;

%Add the gravitational term
    
%     if strcmp(keygravity,'y')
% 
%         RHS_old = RHS_old + Gt;
%         P_new = RHS_old;
%   
%     end
    
%A_old = M_old;
    
% Picard iteration
if strcmp(acel,'FPI')
[p,flowrate,flowresult]=ferncodes_iterpicardNLFVH(M_old,RHS_old,parameter,...
    p_old,nflagface,wells,mobility,weightDMP,gravelem,gravpoint,gravface);
% Iteração de Picard com acelerador de Anderson
elseif strcmp(acel,'AA')
[p,flowrate,flowresult]=ferncodes_iterpicardANLFVH(M_old,RHS_old,...
    parameter,p_old,nflagface,wells,mobility,weightDMP,w,s,nflagno,0,gravelem,gravpoint,gravface);

end
end