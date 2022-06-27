function [p,flowrate,flowresult]=ferncodes_iterpicardANLFVPP2(M_old,RHS_old,...
    parameter,w,s,p_old,nflag,wells,mobility,weightDMP,contnorm)

%% calculo do residuo Inicial
R0=norm(M_old*p_old-RHS_old);
 %% Acelerador de Anderson
   [L,U] = ilu(M_old,struct('type','ilutp','droptol',1e-6));
       
  [p_oldold,fl1,rr1,it1,rv1]=gmres(M_old,RHS_old,10,1e-9,1000,L,U);
  [p,erro,iter]=ferncodes_andersonacc2(p_oldold,1e-6,parameter,w,s,...
    nflag,weightDMP,wells,mobility,R0,contnorm);
%Message to user:
fprintf('\n Iteration number, iterations = %d \n',iter)
fprintf('\n Residual error, error = %d \n',erro)
disp('>> The Pressure field was calculated with success!');
[pinterp]=ferncodes_pressureinterpNLFVPP(p,nflag,w,s);
%Get the flow rate (Diamond)
[flowrate,flowresult]=ferncodes_flowrateNLFVPP(p, pinterp, parameter,mobility);

%Message to user:
disp('>> The Flow Rate field was calculated with success!');

end