function [grav,gravelem,gravpoint,gravface]=PLUG_Gfunction
%function [grav]=PLUG_Gfunction
global  numcase g elem centelem bedge inedge coord keygravity cte

grav=zeros(size(elem,1),3);
gravelem=zeros(size(elem,1),1);
gravpoint=zeros(size(coord,1),1);
gravface=zeros(size(bedge,1)+size(inedge,1),1);
switch numcase
    
    case 18
        for i=1:size(elem,1)
            
            grav(i,:)=g';
            
        end
    
        
        
        
    % 4.2 Numerical Examples: incompressible flow in a unit square domain.
    % The domain has a discontinuity line and the gravitational forces are 
    % given as a linear combination of two contributions
    % STARNONI et al. Consistent MPFA Discretization for Flow in the
    % Presence of Gravity. 2019.
    
    % 4.2 Test 1 where a1=1 e a2=0
    case 421
        
        % Para r=0;
        delta=0.5;
        
        h1=10;
        h2=1;
        
        for i=1: size(elem,1)
            
            y = centelem(i,2);
            
                       
            if(y>delta)
                % gravity vector 
                grav(i,1)=0;
                grav(i,2)=-h1;
                
                % gravity through the edge element
                gravelem(i)= h1*y - cte;
                
            else
                % gravity vector 
                grav(i,1)=0;
                grav(i,2)=-h2;
                
                % gravity through the edge element
                gravelem(i)= h2*y - cte;
            end
                
        end
            
        for j=1: size(coord, 1)

            %Define "x" and "y"
            y2 = coord(j,2);

            if(y2>delta)
                gravpoint(j)= -h1*y2 - cte;
            else
                gravpoint(j)= -h2*y2 - cte;
            end

        end
        
        for j=1:size(bedge,1)+size(inedge,1)
            %Define "x" and "y"
            if j<=size(bedge,1)
                v1=bedge(j,1);
                v2=bedge(j,2);

                a=0.5*(coord(v1,:)+coord(v2,:));
                y11=a(1,2);
            else
                v1=inedge(j-size(bedge,1),1);
                v2=inedge(j-size(bedge,1),2);

                a=0.5*(coord(v1,:)+coord(v2,:));
                y11=a(1,2);
            end

            if(y11>delta)
                gravface(j)= -h1*y11 - cte;
            else
                gravface(j)= -h2*y11 - cte;
            end

        end
             
    
        
        
    % 4.2 Test 2 where a1=0 e a2=1    
    case 422
        
        for i=1: size(elem,1)
            
            x = centelem(i,1);
            y = centelem(i,2);
            
            % gravity vector 
            grav(i,1)= cos(x)*cos(y);
            grav(i,2)=-sin(x)*sin(y);
            
            % gravity through the edge element
            
            gravelem(i,1)=-sin(x)*cos(y)-1+cte;

            
        end
        
        
        for j=1:size(coord,1)
            %Define "x" and "y"
            x1 = coord(j,1);
            y1 = coord(j,2);
            % parametro segundo  Starnoni

            % solucao analitica
            gravpoint(j,1)= -sin(x1)*cos(y1)-1+cte;
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
            gravface(j,1)= -sin(x2)*cos(y2)-1+cte;

        end
            
        
    % 4.2 Test 3 where a1=1 e a2=1
    case 423
        
        % Para r=0;
        delta=0.5;
        
        h1=1;
        h2=10;
        
        for i=1: size(elem,1)
            
            x = centelem(i,1);
            y = centelem(i,2);
                                   
            if(y>delta)
                % gravity vector 
                grav(i,1)= cos(x)*cos(y);
                grav(i,2)=-h1-sin(x)*sin(y);
                
                % gravity through the edge element
                gravelem(i)= -11 + h1*y - sin(x)*cos(y) - cte;

            else
                % gravity vector 
                grav(i,1)= cos(x)*cos(y);
                grav(i,2)=-h2-sin(x)*sin(y);
                
                % gravity through the edge element
                gravelem(i)= -6.5 + h2*y - sin(x)*cos(y) - cte;

            end
                
        end
        
        for j=1: size(coord, 1)

            %Define "x" and "y"
            x1 = coord(j,1);
            y1 = coord(j,2);
            % parametro segundo  Starnoni

            % solucao analitica
            if(y1>delta)
                gravpoint(j)= -11 + h1*y1 - sin(x1)*cos(y1) - cte;
            else
                gravpoint(j)= -6.5 + h2*y1 - sin(x1)*cos(y1) - cte;
            end

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

            if(y2>delta)
                gravface(j)= -11+h1*y2 - sin(x2)*cos(y2) - cte;
            else
                gravface(j)= -6.5+h2*y2 - sin(x2)*cos(y2) - cte;
            end
        end
        

    % 4.2 Test 4 where a1=1 e a2=100
    case 424
        
        % Para r=0;
        delta=0.5;
        
        h1=1;
        h2=10;
        
        for i=1: size(elem,1)
            
            x0 = centelem(i,1);
            y0 = centelem(i,2);
                                   
            if(y0>delta)
                grav(i,1)= 100*cos(x0)*cos(y0);
                grav(i,2)=-h1-100*sin(x0)*sin(y0);
            else
                grav(i,1)= 100*cos(x0)*cos(y0);
                grav(i,2)=-h2-100*sin(x0)*sin(y0);
            end
                
        end
        
        if strcmp(keygravity,'c')
            
            for i=1: size(elem,1)
            
                %Define "x" and "y"
                x1 = centelem(i,1);
                y1 = centelem(i,2);
                
                if(y1>delta)
                    gravelem(i)= -h1*y1 + 100*sin(x1)*cos(y1) - cte;
                else
                    gravelem(i)= -h2*y1 + 100*sin(x1)*cos(y1) - cte;
                end
    
            end
            
            for j=1: size(coord, 1)
               
                %Define "x" and "y"
                x2 = coord(j,1);
                y2 = coord(j,2);
                % parametro segundo  Starnoni
                
                % solucao analitica
                if(y2>delta)
                    gravpoint(j)= -h1*y1 + 100*sin(x2)*cos(y2) - cte;
                else
                    gravpoint(j)= -h2*y1 + 100*sin(x2)*cos(y2) - cte;
                end
                
            end
        
            for j=1:size(bedge,1)+size(inedge,1)
                %Define "x" and "y"
                if j<=size(bedge,1)
                    v1=bedge(j,1);
                    v2=bedge(j,2);
                    
                    a=0.5*(coord(v1,:)+coord(v2,:));
                    x11=a(1,1);
                    y11=a(1,2);
                else
                    v1=inedge(j-size(bedge,1),1);
                    v2=inedge(j-size(bedge,1),2);
                    
                    a=0.5*(coord(v1,:)+coord(v2,:));
                    x11=a(1,1);
                    y11=a(1,2);
                end
                
                if(y11>delta)
                    gravface(j)= -h1*y11 + 100*sin(x11)*cos(y11) - cte;
                else
                    gravface(j)= -h2*y11 + 100*sin(x11)*cos(y11) - cte;
                end
            end
            
        end
        
end

end