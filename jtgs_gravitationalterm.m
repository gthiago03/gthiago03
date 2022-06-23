function [Gt] = jtgs_gravitationalterm(overedgecoord, kmap)

    global inedge bedge centelem elem g normals

    dj = distancej(overedgecoord);
    
    for iface=1:size(bedge,1)+size(inedge,1)

        K_dwh = kdweightedharmonic(dj(iface,:), kmap);
        g_wa = gweightarithmetic(dj(iface,:));
        normal = [normals(iface,1) normals(iface,2)];
        
        g_edge(iface,1)=dot(normal,K_dwh*g_wa);
        g_edge(iface,2)=dj(iface,3);
        g_edge(iface,3)=dj(iface,4);            

    end 

    
    Gt = zeros(size(elem,1),1);
    
    for el=1:size(elem,1)
        
        ind = find(g_edge(:,2)==el | g_edge(:,3)==el);
        

        
        for ij=1:size(ind,1)
            
            Gt(el) = Gt(el) + g_edge(ind(ij));
            
        end
        
    end        
        
    
    %--------------------------------------------------------------------------
    %FUNCTION "distancej"
    %--------------------------------------------------------------------------

    %This function calculates the distance between the center of the analyzed 
    %edge and the center of the elements on its left and right. 
    
    function [dj] = distancej(overedgecoord)
                    
        for iface=1:size(bedge,1)
            
            lel=bedge(iface,3);
            djlef=overedgecoord(iface,:)-centelem(lel,:);

            dj(iface,1)=norm(djlef);
            dj(iface,2)=0;
            dj(iface,3)=lel;
            dj(iface,4)=0;
            
        end
        
        for iface=1:size(inedge,1)

            lel=inedge(iface,3); %elemento da esquerda
            rel=inedge(iface,4); %elemento da direita

            djlef=overedgecoord(size(bedge,1)+iface,:)-centelem(lel,:);
            djrig=overedgecoord(size(bedge,1)+iface,:)-centelem(rel,:);

            dj(iface+size(bedge,1),1)=norm(djlef);
            dj(iface+size(bedge,1),2)=norm(djrig);
            dj(iface+size(bedge,1),3)=lel;
            dj(iface+size(bedge,1),4)=rel;

        end
    
    end
    
    %--------------------------------------------------------------------------
    %FUNCTION "dweightdharmonic"
    %--------------------------------------------------------------------------

    %This function calculates the d-weighted harmonic average of the 
    %permeability tensors between the left and right elements 
    
    function [K_dwh] = kdweightedharmonic(dj, kmap)
        
        if(dj(4)==0)
        
            lel=dj(3);
            
            Klel=zeros(2,2); % left element tensor
            
            % left element tensor
            
            Klel(1,1)=kmap(elem(lel,5),2);
            Klel(1,2)=kmap(elem(lel,5),3);
            Klel(2,1)=kmap(elem(lel,5),4);
            Klel(2,2)=kmap(elem(lel,5),5);

            K_dwh=inv(dj(1)*inv(Klel));         
            
        else

            lel=dj(3);
            rel=dj(4);

            Klel=zeros(2,2); % left element tensor
            Krel=zeros(2,2); % right element tensor

            % left element tensor

            Klel(1,1)=kmap(elem(lel,5),2);
            Klel(1,2)=kmap(elem(lel,5),3);
            Klel(2,1)=kmap(elem(lel,5),4);
            Klel(2,2)=kmap(elem(lel,5),5);

            % right element tensor

            Krel(1,1)=kmap(elem(rel,5),2);
            Krel(1,2)=kmap(elem(rel,5),3);
            Krel(2,1)=kmap(elem(rel,5),4);
            Krel(2,2)=kmap(elem(rel,5),5);
        
            
            K_dwh=inv(dj(1)*inv(Klel)+dj(2)*inv(Krel));
            
        
        end
    end

    %--------------------------------------------------------------------------
    %FUNCTION "gweightedharmonic"
    %--------------------------------------------------------------------------

    %This function calculates the weighted arithmetic average of the edge 
    %gravity vectors 
    function [g_wa] = gweightarithmetic(dj)
        
        if(dj(4)==0)
            
            lel=dj(3);
            
            g_wa=dj(1)*[g(1) g(2)]';          
            
        else 

            lel=inedge(3);
            rel=inedge(4);

            g_wa=dj(1)*[g(1) g(2)]'+dj(2)*[g(1) g(2)]';

            
        end        
        
    end

end