%It is called by the function "ferncodes_Pre_LPEW_1.m"

function [Kt1,Kt2,Kn1,Kn2] = ferncodes_Ks_Interp_LPEW1(O,T,Qo,kmap,No,...
    mobility,Sw,V)
%Retorna os K(n ou t) necess�rios para a obten��o dos weights. kmap � a
%matriz de permeabilidade; Ex: Kt1->linhaN=Kt1(cellN);

global bedge inedge esurn2 esurn1 phasekey visc elem numcase

%Initialize "bedgesize" and "inedgesize"
bedgesize = size(bedge,1);
inedgesize = size(inedge,1);

%Verifica quantos s�o os elementos em torno do n� "ni".%
N_element_No=esurn2(No+1)-esurn2(No);

%Prealoca��o das matrizes.%

Kt1=zeros(N_element_No,2); %As duas colunas correspondem a i=1 e i=2.
Kt2=zeros(N_element_No,2);
Kn1=zeros(N_element_No,2);
Kn2=zeros(N_element_No,2);
K=zeros(3);
K1=zeros(3);

R=[0 1 0; -1 0 0;0 0 0];

%--------------------------------------------------------------------------
%Verify if it is an One-Phase or Two-Phase case (FACE "mobility" definit.)
    
%One-Phase case:
if phasekey == 1
    mobility = ones(bedgesize + inedgesize,1);
%Two-Phase case
else
    %It sums the water and oil mobilities (came from Marcio's code)
    mobility = sum(mobility,2);
    %Reorganize the position of total mobility for it be used in Fernando's 
    %code
    %Put in "auxmobility" the mobilities calculated in the boundary edges.
    auxmobility = mobility(1:bedgesize);
    %Put the mobilities calculated in the internal edges in first place
    mobility(1:inedgesize) = mobility(bedgesize + 1:length(mobility));
    %Complete the mobily filling
    mobility(inedgesize + 1:length(mobility)) = auxmobility; 
end  %End of IF

%--------------------------------------------------------------------------

%constru��o do tensor permeabilidade.%

for k=1:N_element_No
    
    j=esurn1(esurn2(No)+k);
    for icont=1:2
        if (size(T,1)==size(O,1))&&(k==N_element_No)&&(icont==2)
                        
            K(1,1)=mobility(V(icont,k,No))*kmap(elem(j,5),2);
            K(1,2)=mobility(V(icont,k,No))*kmap(elem(j,5),3);
            K(2,1)=mobility(V(icont,k,No))*kmap(elem(j,5),4);
            K(2,2)=mobility(V(icont,k,No))*kmap(elem(j,5),5);
            
            Kn1(k,icont)=((R*(T(1,:)-Qo)')'*K*(R*(T(1,:)-Qo)'))/(norm(T(1,:)-Qo)^2);
            Kt1(k,icont)=((R*(T(1,:)-Qo)')'*K*((T(1,:)-Qo)'))/(norm(T(1,:)-Qo)^2);
        else
            K(1,1)=mobility(V(icont,k,No))*kmap(elem(j,5),2);
            K(1,2)=mobility(V(icont,k,No))*kmap(elem(j,5),3);
            K(2,1)=mobility(V(icont,k,No))*kmap(elem(j,5),4);
            K(2,2)=mobility(V(icont,k,No))*kmap(elem(j,5),5);
                     
            Kn1(k,icont)=((R*(T(k+icont-1,:)-Qo)')'*K*(R*(T(k+icont-1,:)-Qo)'))/(norm(T(k+icont-1,:)-Qo)^2);
            Kt1(k,icont)=((R*(T(k+icont-1,:)-Qo)')'*K*((T(k+icont-1,:)-Qo)'))/(norm(T(k+icont-1,:)-Qo)^2);
        end
    end
    
    %----------------------------------------------------------------------
    %Verify if it is an One-Phase or Two-Phase case (ELEMENT "mobility")
 
    %Verify if it is an One-Phase or Two-Phase case
    %One-Phase case:
    if phasekey == 1
        L22 = 1;
    %Two-Phase case
    else
        %Call "twophasevar.m"
        [null1,null2,null3,krw,kro,] = twophasevar(Sw(j),numcase);
 
        L22 = krw/visc(1) + kro/visc(2);   
    end  %End of IF
    
    %------------------------- Tensores ----------------------------------%
    K1(1,1)=L22*kmap(elem(j,5),2);
    K1(1,2)=L22*kmap(elem(j,5),3);
    K1(2,1)=L22*kmap(elem(j,5),4);
    K1(2,2)=L22*kmap(elem(j,5),5);
    
    % calculo dos outros K(n ou t) no paper � denotado com um "~" na parte
    % inferior de K
    for icont=0:1
        if (size(T,1)==size(O,1))&&(k==N_element_No)&&(icont==1)
            
            Kn2(k,icont+1)=((R*(O(k,:)-T(1,:))')'*K1*(R*(O(k,:)-T(1,:))'))/(norm(O(k,:)-T(1,:))^2);
            Kt2(k,icont+1)=((R*(O(k,:)-T(1,:))')'*K1*((O(k,:)-T(1,:))'))/(norm(O(k,:)-T(1,:))^2);
            
        else
            
            Kn2(k,icont+1)=((R*(O(k,:)-T(k+icont,:))')'*K1*(R*(O(k,:)-T(k+icont,:))'))/(norm(O(k,:)-T(k+icont,:))^2);
            Kt2(k,icont+1)=((R*(O(k,:)-T(k+icont,:))')'*K1*((O(k,:)-T(k+icont,:))'))/(norm(O(k,:)-T(k+icont,:))^2);
        end
    end
    
end

end

