function [grav,gravelem,gravface]=PLUG_Gfunction
global keygravity numcase g elem centelem bedge inedge coord

grav=zeros(size(elem,1),3);

switch numcase
    
    case 18
        for i=1:size(elem,1)
            
            grav(i,:)=g';
            
        end
    
        
%         for i = 1:size(centelem,1)
%             %Define "x" and "y"
%             x = centelem(i,1);
%             y = centelem(i,2);
%             % solucao analitica foi calculado usando pag. 385
%             % Calculo II Tom Apostol
%             u(i,1)= -sind(x)*cosd(y);
%             gravelem(i,1)=sind(x)*cosd(y);
%             % gravidade
%             grav(i,:)=[cosd(x)*cosd(y) -sind(x)*sind(y)];
%             
%         end
        
        
    % 4.2 Numerical Examples: incompressible flow in a unit square domain.
    % The domain has a discontinuity line and the gravitational forces are 
    % given as a linear combination of two contributions
    % STARNONI et al. Consistent MPFA Discretization for Flow in the
    % Presence of Gravity. 2019.
    
    % 4.2 Test 1 where a1=1 e a2=0
    case 421
        
        r=0;
        s=1-r;
        delta=0.5;
        F=r*centelem(1)+s*centelem(2);
            
        if(F>delta)
            h=1;
        else
            h=10;
        end

        for i=1: size(elem,1)
            
            grav(i,1)= -(r*h)/sqrt(r^2+s^2);
            grav(i,2)= -(s*h)/sqrt(r^2+s^2);
        
        end
        
        
    % 4.2 Test 2 where a1=0 e a2=1    
    case 422
        
        for i=1: size(elem,1)
            
            grav(i,1)= cosd(centelem(i,1))*cosd(centelem(i,2));
            grav(i,2)= -sind(centelem(i,1))*sind(centelem(i,2));
            
                %Define "x" and "y"
                x = centelem(i,1);
                y = centelem(i,2);
                % solucao analitica foi calculado usando pag. 385
                % Calculo II Tom Apostol
               
                gravelem(i,1)=sind(x)*cosd(y);
           
 
        end
        for j=1:size(bedge,1)+size(inedge,1)
            %Define "x" and "y"
            if j<=size(bedge,1)
                v1=bedge(j,1);
                v2=bedge(j,2);
                
                a=0.5*(coord(v1,:)+coord(v2,:));
                x2=a(1,1);
                y2=a(1,2);
            else
                v1=inedge(j-size(bedge,1),1);
                v2=inedge(j-size(bedge,1),2);
                
                a=0.5*(coord(v1,:)+coord(v2,:));
                x2=a(1,1);
                y2=a(1,2);
            end
            gravface(j,1)= sind(x2)*cosd(y2);
        end
        
    % 4.2 Test 3 where a1=1 e a2=1
    case 423
        
        r=0;
        s=1-r;
        delta=0.5;
        F=r*centelem(1)+s*centelem(2);
            
        if(F>delta)
            h=1;
        else
            h=10;
        end

        for i=1: size(elem,1)
            
            grav(i,1)= -(r*h)/sqrt(r^2+s^2)+cosd(centelem(i,1))*cosd(centelem(i,2));
            grav(i,2)= -(s*h)/sqrt(r^2+s^2)-sind(centelem(i,1))*sind(centelem(i,2));
        
        end
        
    % 4.2 Test 4 where a1=1 e a2=100
    case 424
        
        r=0;
        s=1-r;
        delta=0.5;
        F=r*centelem(1)+s*centelem(2);
            
        if(F>delta)
            h=1;
        else
            h=10;
        end

        for i=1: size(elem,1)
            
            grav(i,1)= -(r*h)/sqrt(r^2+s^2)+100*cosd(centelem(i,1))*cosd(centelem(i,2));
            grav(i,2)= -(s*h)/sqrt(r^2+s^2)-100*sind(centelem(i,1))*sind(centelem(i,2));
        
        end
        
end

end