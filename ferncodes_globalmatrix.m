%It is called by "ferncodes_solvepressure.m"

function [M,I] = ferncodes_globalmatrix(w,s,Kde,Ded,Kn,Kt,Hesq,...
    nflag)
%Define global variables:
global coord elem esurn1 esurn2 bedge inedge centelem bcflag ...
    phasekey;

%-----------------------inicio da rOtina ----------------------------------%
%Constr�i a matriz global.

%Initialize "bedgesize" and "inedgesize"
bedgesize = size(bedge,1);
inedgesize = size(inedge,1);

%Initialize "M" (global matrix) and "I" (known vector)
M = sparse(size(elem,1),size(elem,1)); %Prealoca��o de M.
I = zeros(size(elem,1),1);

% % contribui��o dos po�os
% if max(wells)~=0
%     sumvol=0;
%     for iw = 1:size(wells,1)
%         
%         if wells(iw,2)==1            % injetor
%             I(wells(iw,1))= 1*elemarea(wells(iw,1));        % injeta um m3 de agua por dia (d)
%             sumvol = sumvol + elemarea(wells(iw,1));
%         end
%     end
%     I = I./sumvol;
% else
    
%Evaluate "bedge"
for ifacont = 1:bedgesize
            
    %Get element on the left:
    lef = bedge(ifacont,3);
    %Get another parameters:
    v0 = coord(bedge(ifacont,2),:) - coord(bedge(ifacont,1),:); %fase.
    v1 = centelem(bedge(ifacont,3),:) - coord(bedge(ifacont,1),:);
    v2 = centelem(bedge(ifacont,3),:) - coord(bedge(ifacont,2),:);
    
    %Calculate the nom of the edge
    nor = norm(coord(bedge(ifacont,1),:) - coord(bedge(ifacont,2),:));

    % Tratamento do n� nos v�rtices 2 e 4%
        
    %Dirichlet Boundary
    if bedge(ifacont,5) < 200
        c1 = nflag(bedge(ifacont,1),2);
        c2 = nflag(bedge(ifacont,2),2);
        
        A = -Kn(ifacont)/(Hesq(ifacont)*norm(v0));
            
        %Preenchimento
        
        M(bedge(ifacont,3),bedge(ifacont,3)) = M(bedge(ifacont,3),...
            bedge(ifacont,3)) - A*(norm(v0)^2);
            
        %!!!!!!!!!!!!!!!!!!!!Vericar no monof�sico o sinal
        I(bedge(ifacont,3)) = I(bedge(ifacont,3)) - ...
            (dot(v2,-v0)*c1 + dot(v1,v0)*c2)*A + ...
            (c2 - c1)*Kt(ifacont); 
            
    %Neumann boundary
    else
        x = logical(bcflag(:,1) == bedge(ifacont,5));
        I(lef) = I(lef) + nor*bcflag(x,2);
    end  %End of IF
end  %End of FOR

% end  %End of IF

% contribui��o nas faces internas
for iface = 1:inedgesize
    
    % press�o prescrita no elemento do po�o injetor
    %Contabiliza as contribui��es do fluxo numa aresta para os elementos %
    %a direita e a esquerda dela.                                        %
    
    M(inedge(iface,3),inedge(iface,3)) = ...
        M(inedge(iface,3),inedge(iface,3)) - Kde(iface);
    M(inedge(iface,3),inedge(iface,4)) = ...
        M(inedge(iface,3),inedge(iface,4)) + Kde(iface);
    M(inedge(iface,4),inedge(iface,4)) = ...
        M(inedge(iface,4), inedge(iface,4)) - Kde(iface);
    M(inedge(iface,4),inedge(iface,3)) = ...
        M(inedge(iface,4), inedge(iface,3)) + Kde(iface);
    
    if nflag(inedge(iface,1),1) < 200
        I(inedge(iface,3)) = I(inedge(iface,3)) - Kde(iface)*...
            Ded(iface)*nflag(inedge(iface,1),2);
        I(inedge(iface,4)) = I(inedge(iface,4)) + Kde(iface)*...
            Ded(iface)*nflag(inedge(iface,1),2);
    end
    if nflag(inedge(iface,2),1) < 200
        I(inedge(iface,3)) = I(inedge(iface,3)) +Kde(iface)*...
            Ded(iface)*nflag(inedge(iface,2),2);
        I(inedge(iface,4)) = I(inedge(iface,4)) - Kde(iface)*...
            Ded(iface)*nflag(inedge(iface,2),2);
    end
    % quando o n� pertece ao contorno de Neumann
    if nflag(inedge(iface,1),1) == 202
        
        I(inedge(iface,3)) = I(inedge(iface,3)) - Kde(iface)*...
            Ded(iface)*s(inedge(iface,1)); %ok
        
        I(inedge(iface,4)) = I(inedge(iface,4)) + Kde(iface)*...
            Ded(iface)*s(inedge(iface,1)); %ok
    end
    if nflag(inedge(iface,2),1) == 202
        
        I(inedge(iface,3)) = I(inedge(iface,3)) + Kde(iface)*...
            Ded(iface)*s(inedge(iface,2)); %ok
        
        I(inedge(iface,4)) = I(inedge(iface,4)) - Kde(iface)*...
            Ded(iface)*s(inedge(iface,2)); %ok
        
    end
    
    %Contabiliza��o das contribui��es dos n�s que n�o est�o na
    %fronteiras de Dirichlet.
    
    if nflag(inedge(iface,1),1) > 200
        for j = 1:(esurn2(inedge(iface,1) + 1) - esurn2(inedge(iface,1)))
            
            post_cont = esurn2(inedge(iface,1)) + j;
            
            M(inedge(iface,3),esurn1(post_cont)) = M(inedge(iface,3),...
                esurn1(post_cont)) + Kde(iface)*Ded(iface)*w(post_cont);
            
            M(inedge(iface,4),esurn1(post_cont)) = M(inedge(iface,4),...
                esurn1(post_cont)) - Kde(iface)*Ded(iface)*w(post_cont);
            
        end
    end
    if nflag(inedge(iface,2),1) > 200
        for j = 1:(esurn2(inedge(iface,2) + 1) - esurn2(inedge(iface,2)))
            
            post_cont = esurn2(inedge(iface,2)) + j;
            
            M(inedge(iface,3), esurn1(post_cont)) = M(inedge(iface,3),...
                esurn1(post_cont)) - Kde(iface)*Ded(iface)*w(post_cont);
            
            M(inedge(iface,4), esurn1(post_cont)) = M(inedge(iface,4),...
                esurn1(post_cont)) + Kde(iface)*Ded(iface)*w(post_cont);
        end
    end
   
end  %End of FOR ("inedge")    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
%     %----------------------------------------------------------------------
%     %Boundary treatment (faces)
%     
%     %The face is inside the domain, however one or both the node(s) may be 
%     %over the boundary.
%     
%     %-----------------------------------------
%     %Verify if some vertex is on the boundary:
%     
%     %Get "vertices"
%     vertices = inedge(iface,1:2);
% 
%     %Se os n�s das arestas estiverem em fronteiras de Dirichlet, suas
%     %contribui��es ser�o contabilizadas logo abaixo.
%     
%     %--------------------------------------
%     %The vertex 1 has a Dirichlet boundary:
% 
%     %It points if the vertex belong to boundary: 
%     pointvtxonbound = logical(vertices(1) == bedge(:,1));
% 
%     %Verify if the vertex is on the boundary
%     if any(pointvtxonbound)
%         %Get the row in "bedge" where this occurs
%         bedgrow = bedgeamount(pointvtxonbound);
% 
%         %There exists a vertex on the boundary and it is a 
%         %Dirichlet boundary:
%         if bedge(bedgrow(1),4) < 200
%             %Define "flagpointer"
%             flagpointer = logical(bcflag(:,1) == bedge(bedgrow(1),4));
%             
%             %################################################
%             %Adapted to generalized Marcio's code
%             knownval = PLUG_bcfunction(vertices(1),flagpointer);  %nflag(bedge(ifacont,1),2);
%             %################################################
% % vertices
% % knownval
%             I(inedge(iface,3)) = I(inedge(iface,3)) - mobonface*Kde(iface)*...
%                 Ded(iface)*knownval;  %It was modified (M�rcio)
%             I(inedge(iface,4)) = I(inedge(iface,4)) + mobonface*Kde(iface)*...
%                 Ded(iface)*knownval;  %It was modified (M�rcio)
%         end  %End of IF
% 
%         %---------------------------------------------
%         %The vertex 1 has a Neumann boundary:
% 
%         % quando o n� pertece ao contorno de Neumann
%         if bedge(bedgrow(1),4) > 200
%             %Define "flagpointer"
%             flagpointer = logical(bcflag(:,1) == bedge(bedgrow(1),4));
%             %Verify if the boundary value is null
%             %NON-Null value (in Fernando's code it was a 202 flag)
%             if bcflag(flagpointer,2) ~= 0
%                 I(inedge(iface,3)) = I(inedge(iface,3)) - ...
%                     mobonface*Kde(iface)*Ded(iface)*s(inedge(iface,1)); %ok
%         
%                 I(inedge(iface,4)) = I(inedge(iface,4)) + ...
%                     mobonface*Kde(iface)*Ded(iface)*s(inedge(iface,1)); %ok
%             %NULL value (== 0)
%             else
%                 for j = 1:(esurn2(inedge(iface,1) + 1) - ...
%                         esurn2(inedge(iface,1)))
%                     post_cont = esurn2(inedge(iface,1)) + j;
%             
%                     M(inedge(iface,3),esurn1(post_cont)) = ...
%                         M(inedge(iface,3),esurn1(post_cont)) + ...
%                         mobonface*Kde(iface)*Ded(iface)*w(post_cont);
%             
%                     M(inedge(iface,4),esurn1(post_cont)) = ...
%                         M(inedge(iface,4),esurn1(post_cont)) - ...
%                         mobonface*Kde(iface)*Ded(iface)*w(post_cont);
%                 end  %End of FOR
%             end  %End of IF
%         end  %End of IF
%     end  %End of IF (the vertex belongs to boundary)
%     
%     %--------------------------------------
%     %The vertex 2 has a Dirichlet boundary:
% 
%     %It points if the vertex belong to boundary: 
%     pointvtxonbound = logical(vertices(2) == bedge(:,1));
%     %The vertex belongs to boundary
%     if any(pointvtxonbound)
%         %Get the row in "bedge" where this occurs
%         bedgrow = bedgeamount(pointvtxonbound);
%     
%         %There exists a vertex on the boundary and it is a Dirichlet boundary:
%         if bedge(bedgrow(1),4) < 200
%             %Define "flagpointer"
%             flagpointer = logical(bcflag(:,1) == bedge(bedgrow(1),4));
%             
%             %################################################
%             %Adapted to generalized Marcio's code
%             knownval = PLUG_bcfunction(vertices(2),flagpointer);  %nflag(bedge(ifacont,1),2);
%             %################################################
% % knownval
% % pause
%             I(inedge(iface,3)) = I(inedge(iface,3)) + mobonface*Kde(iface)*...
%                 Ded(iface)*knownval;
%             I(inedge(iface,4)) = I(inedge(iface,4)) - mobonface*Kde(iface)*...
%                 Ded(iface)*knownval;
%         end  %End of IF 
%     
%         %------------------------------------
%         %The vertex 2 has a Neumann boundary:
%     
%         % quando o n� pertece ao contorno de Neumann
%         if bedge(bedgrow(1),4) > 200
%             %Define "flagpointer"
%             flagpointer = logical(bcflag(:,1) == bedge(bedgrow(1),4));
%             %Verify if the boundary value is null
%             %NON-Null value (in Fernando's code it was a 202 flag)
%             if bcflag(flagpointer,2) > 0
%                 %Left contribution
%                 I(inedge(iface,3)) = I(inedge(iface,3)) + ...
%                     mobonface*Kde(iface)*Ded(iface)*s(inedge(iface,2)); %ok
%                 %Right contribution
%                 I(inedge(iface,4)) = I(inedge(iface,4)) - ...
%                     mobonface*Kde(iface)*Ded(iface)*s(inedge(iface,2)); %ok
%             %NULL value (== 0)
%             else
%                 for j = 1:(esurn2(inedge(iface,2) + 1) - ...
%                         esurn2(inedge(iface,2))),
%                     post_cont = esurn2(inedge(iface,2)) + j;
%             
%                     M(inedge(iface,3),esurn1(post_cont)) = ...
%                         M(inedge(iface,3),esurn1(post_cont)) - ...
%                         mobonface*Kde(iface)*Ded(iface)*w(post_cont);
%             
%                     M(inedge(iface,4),esurn1(post_cont)) = ...
%                         M(inedge(iface,4),esurn1(post_cont)) + ...
%                         mobonface*Kde(iface)*Ded(iface)*w(post_cont);
%                 end  %End of FOR
%             end  %End of IF
%         end  %End of IF
%     end  %End of IF (the vertex 2 belongs to boundary)
%  end  %End of FOR ("inedge" faces)


% 
% % adequa��o da matriz nos po�os produtores
% if max(wells)~=0
%     for iw = 1:size(wells,1)
%         if wells(iw,2) == 2 %produtor
%             M(wells(iw,1),:)=0*M(wells(iw,1),:);
%             M(wells(iw,1),wells(iw,1))=1;
%             I(wells(iw,1)) = 0;
%         end  %End of IF
%     end  %End of FOR
% end  %End of IF


