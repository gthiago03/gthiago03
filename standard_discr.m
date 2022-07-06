function [nflagface]=gravitation(kmap,vec_gravelem, vec_gravface, vec_gravpoint, nflagface)
%atribui pressões aos pontos e faces caso o método de discretização do
%termo gravitacional seja o padrão

    global bedge elem centelem centface coord normals
    
    Klef=zeros(2,2);
    
    for ifacont=1 : size(bedge,1)

        lef = bedge(ifacont,3);
        
        % tensor de permeabilidade do elemento a esquerda
        Klef(1,1)=kmap(elem(lef,5),2);
        Klef(1,2)=kmap(elem(lef,5),3);
        Klef(2,1)=kmap(elem(lef,5),4);
        Klef(2,2)=kmap(elem(lef,5),5);
        
        if(bedge(ifacont,5) > 200)

%             kgradpx_p = Klef(1,1)*vec_gravpoint(ifacont,1)+Klef(1,2)*vec_gravpoint(ifacont,2);
%             kgradpy_p = Klef(2,1)*vec_gravpoint(ifacont,1)+Klef(2,2)*vec_gravpoint(ifacont,2);
%             
%             kgx_p = Klef(1,1)*vec_gravpoint(ifacont,1)+Klef(1,2)*vec_gravpoint(ifacont,2);
%             kgy_p = Klef(2,1)*vec_gravpoint(ifacont,1)+Klef(2,2)*vec_gravpoint(ifacont,2);
            
            
            kgradpx = Klef(1,1)*vec_gravface(ifacont,1)+Klef(1,2)*vec_gravface(ifacont,2);
            kgradpy = Klef(2,1)*vec_gravface(ifacont,1)+Klef(2,2)*vec_gravface(ifacont,2);
            
            kgx = Klef(1,1)*vec_gravface(ifacont,1)+Klef(1,2)*vec_gravface(ifacont,2);
            kgy = Klef(2,1)*vec_gravface(ifacont,1)+Klef(2,2)*vec_gravface(ifacont,2);
          
            
            flux_ex_du = - (abs(normals(ifacont,1))*kgradpx + abs(normals(ifacont,2))*kgradpy);
            flux_ex_g = abs(normals(ifacont,1))*kgx + abs(normals(ifacont,2))*kgy;
            flux_ex = flux_ex_du + flux_ex_g;
            
            if(centface(ifacont,2)==0)
                nflagface(ifacont,2) = (flux_ex - flux_ex_g)*(-1);
            else
                nflagface(ifacont,2) = flux_ex - flux_ex_g;
            end
            
           % if(coord(ifacont,2)==0)
           %     nflag(ifacont,2) = (flux_ex - flux_ex_g)*(-1);
            
        end
        
    end
    
end

% nflag(5,2) = (nflagface(1,2)+nflagface(2,2))/2;
% nflag(6,2) = (nflagface(2,2)+nflagface(3,2))/2;
% nflag(7,2) = (nflagface(3,2)+nflagface(4,2))/2;
% nflag(11,2) = (nflagface(9,2)+nflagface(10,2))/2;
% nflag(12,2) = (nflagface(10,2)+nflagface(11,2))/2;
% nflag(13,2) = (nflagface(11,2)+nflagface(12,2))/2;