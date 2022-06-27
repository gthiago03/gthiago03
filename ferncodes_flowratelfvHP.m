function [flowrate,flowresult]=ferncodes_flowratelfvHP(parameter,weightDMP,pinterp,p)
global inedge coord bedge bcflag centelem 
%Initialize "bedgesize" and "inedgesize"
bedgesize = size(bedge,1);
inedgesize = size(inedge,1);

%Initialize "flowrate" and "flowresult"
flowrate = zeros(bedgesize + inedgesize,1);
flowresult = zeros(size(centelem,1),1);
    
for ifacont=1:bedgesize
        
    lef=bedge(ifacont,3);
    
    normcont=norm(coord(bedge(ifacont,1),:)-coord(bedge(ifacont,2),:));
    
    if bedge(ifacont,5)>200
        x=bcflag(:,1)==bedge(ifacont,5);
        r=find(x==1);
        %flowrate(ifacont,1)=normcont*bcflag(r,2);% testes feitos em todos os
        %problemas monof�sico
        flowrate(ifacont,1)=-normcont*bcflag(r,2);% problema de buckley leverett Bastian
    else
        facelef1=parameter(1,3,ifacont);
        facelef2=parameter(1,4,ifacont);
        
        flowrate(ifacont,1)= normcont*((parameter(1,1,ifacont)+parameter(1,2,ifacont))*p(lef)-...
            parameter(1,1,ifacont)*pinterp(facelef1)-parameter(1,2,ifacont)*pinterp(facelef2));
    end
    %Attribute the flow rate to "flowresult"
    %On the left:
    flowresult(lef) = flowresult(lef) + flowrate(ifacont);
    
end
% Montagem da matriz global

for iface=1:inedgesize
    
    lef=inedge(iface,3);
    rel=inedge(iface,4);
    %Determina��o dos centr�ides dos elementos � direita e � esquerda.%
    vd1=coord(inedge(iface,2),:)-coord(inedge(iface,1),:);
    norma=norm(vd1);
    ifactual=iface+size(bedge,1);
    % Calculo das contribui��es do elemento a esquerda
    mulef=weightDMP(ifactual-size(bedge,1),1);
    murel=weightDMP(ifactual-size(bedge,1),2);
   
    % os n�s que conforman os pontos de interpola��o no elemento a esquerda
    auxfacelef1=parameter(1,3,ifactual);
    auxfacelef2=parameter(1,4,ifactual);
    % os n�s que conforman os pontos de interpola��o no elemento a direita
    auxfacerel1=parameter(2,3,ifactual);
    auxfacerel2=parameter(2,4,ifactual);
    % calculo dos fluxo parcial a esquerda
    fluxesq=norma*((parameter(1,1,ifactual)+parameter(1,2,ifactual))*p(lef)-...
                   parameter(1,1,ifactual)*pinterp(auxfacelef1)-parameter(1,2,ifactual)*pinterp(auxfacelef2));
    % calculo dos fluxo parcial a direita        
    fluxdireit=norma*((parameter(2,1,ifactual)+parameter(2,2,ifactual))*p(rel)-...
                       parameter(2,1,ifactual)*pinterp(auxfacerel1)-parameter(2,2,ifactual)*pinterp(auxfacerel2));
    % calculo do fluxo unico na face
    flowrate(iface+size(bedge,1),1)=(murel*fluxesq-mulef*fluxdireit);
    
    %Attribute the flow rate to "flowresult"
    %On the left:
    flowresult(lef) = flowresult(lef) + flowrate(bedgesize + iface);  
    %On the right:
    flowresult(rel) = flowresult(rel) - flowrate(bedgesize + iface);
    
end
end