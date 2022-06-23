function [ep, eq, pi_exact] = convergence(pi)

    global numcase elem centelem elemarea u0
    pi_exact = zeros(size(elem,1),1);
    
    switch numcase
    
        case 421
            
            % Para r=0;
            delta=0.5;
        
            h1=10;
            h2=1;
                    
            for i=1: size(elem,1)
                
                y=centelem(i,2);
                
                if(y>delta)
                    pi_exact(i)=11-h1*y + cte;
                else
                    pi_exact(i)=6.5-h2*y + cte;
                end
                
            end
 
        
        case 422
            
           
            % Calculates exact pressure field
            for i=1: size(elem,1)
           
                x = centelem(i,1);
                y = centelem(i,2);
                
                pi_exact(i,1)= sin(x)*cos(y) + u0;
                
            end
        
        case 423
            
            % Calculates exact pressure field
            % Para r=0;
            delta=0.5;
        
            h1=10;
            h2=1;
        
            for i=1: size(elem,1)
                
                x=centelem(i,1);
                y=centelem(i,2);
                
                if(y>delta)
                    pi_exact(i)=11-h1*y + sin(x)*cos(y) + cte;
                else
                    pi_exact(i)=6.5-h2*y + sin(x)*cos(y) + cte;
                end
                
            end
            
        case 424
            
            % Calculates exact pressure field
            % Para r=0;
            delta=0.5;
        
            h1=10;
            h2=1;
        
            for i=1: size(elem,1)
                
                x=centelem(i,1);
                y=centelem(i,2);
                
                if(y>delta)
                    pi_exact(i)=-h1*y + 100*sin(x)*cos(y)+11 + cte;
                else
                    pi_exact(i)=-h2*y + 100*sin(x)*cos(y)+6.5 + cte;
                end
                
            end    
            
        otherwise
            
            ep=0;
            eq=0;
            
    end
    
                       
            sum1=0;
            sum2=0;
            
            for i=1: size(elem,1)
                
                sum1 = sum1 + elemarea(i)*(pi(i)-pi_exact(i))^2;
                sum2 = sum2 + elemarea(i)*pi_exact(i)^2;
                
            end
    
    
            %Eq. 42 Starnoni et al.
            ep=sum1^(1/2)/sum2^(1/2);
            
            
            eq=0;



end