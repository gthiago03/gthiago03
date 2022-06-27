%It is called by "ferncodes_solverpressure.m"

function [w,s] = ferncodes_Pre_LPEW_2_con(dmap,N)
%Define global parameters
global coord nsurn1 nsurn2 bcflagc bedge inedge

% Retorna todos os par�metros necess�rios �s express�es dos fluxos.
apw = ones(size(coord,1),1);
r = zeros(size(coord,1),2);

for y = 1:size(coord,1),
    No = y;
    % calculos dos vetores O, P, T, Q
    [O,P,T,Qo] = OPT_Interp_LPEW(No);
    % calculo dos angulos
    [ve2,ve1,theta2,theta1] = angulos_Interp_LPEW2(O,P,T,Qo,No);
    % calculo dos netas
    [neta] = netas_Interp_LPEW(O,P,T,Qo,No);
    
    % calculo dos Ks
    [Kt1,Kt2,Kn1,Kn2] = ferncodes_Ks_Interp_LPEW2(O,T,Qo,dmap,No);
    
    % calculo dos lamdas
    [lambda,r] = Lamdas_Weights_LPEW2(Kt1,Kt2,Kn1,Kn2,theta1,theta2,ve1,...
        ve2,neta,P,O,Qo,No,T,r);
    % calculo dos pesos
    for k = 0:size(O,1) - 1,
        w(apw(No) + k) = lambda(k + 1)/sum(lambda); %Os pesos fazem sentido%%%%%%%%%%%
    end
    
    apw(No + 1) = apw(No) + size(O,1);
    % interpola�ao das press�es nos contornos de Neumann
    vetor = nsurn1(nsurn2(No) + 1:nsurn2(No + 1));
    comp1 = N(No,1);
    comp2 = N(No,length(vetor));
    if comp1 > size(inedge,1) && comp2 > size(inedge,1)
        a = bcflagc(:,1) == bedge(comp1 - size(inedge,1),7);
        s1 = find(a == 1);
        b = bcflagc(:,1) == bedge(comp2 - size(inedge,1),7);
        s2 = find(b == 1);
        
        s(No,1) = -(1/sum(lambda))*(r(No,1)*bcflagc(s1,2) + ...
            r(No,2)*bcflagc(s2,2));
    end  %End of IF
end  %End of FOR




