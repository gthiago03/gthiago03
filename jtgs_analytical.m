function [pi_exact] = jtgs_analytical()

    global numcase elem centelem cte
    pi_exact = zeros(size(elem,1),1);
    
    switch numcase
    
        case 421
            
            % Para r=0;
            delta=0.5;
        
            h1=1;
            h2=10;
        
            for i=1: size(elem,1)
                          
                if(centelem(i,2)>delta)
                    pi_exact(i)=h1*delta + cte;
                else
                    pi_exact(i)=h2*delta + cte;
                end
                
            end
 
        
        case 422
            
            % Calculates exact pressure field
            for i=1: size(elem,1)
           
                pi_exact(i,1)= -sin(centelem(i,1))*cos(centelem(i,2)) + cte;
                
            end
        
        case 423
            
            % Calculates exact pressure field
            % Para r=0;
            delta=0.5;
        
            h1=1;
            h2=10;
        
            for i=1: size(elem,1)
                          
                if(centelem(i,2)>delta)
                    pi_exact(i)=h1*delta -sin(centelem(i,1))*cos(centelem(i,2)) + cte;
                else
                    pi_exact(i)=h2*delta -sin(centelem(i,1))*cos(centelem(i,2)) + cte;
                end
                
            end
            
        case 424
            
            % Calculates exact pressure field
            % Para r=0;
            delta=0.5;
        
            h1=1;
            h2=10;
        
            for i=1: size(elem,1)
                          
                if(centelem(i,2)>delta)
                    pi_exact(i)=h1*delta -100*sin(centelem(i,1))*cos(centelem(i,2)) + cte;
                else
                    pi_exact(i)=h2*delta -100*sin(centelem(i,1))*cos(centelem(i,2)) + cte;
                end
                
            end    
            
    end
end