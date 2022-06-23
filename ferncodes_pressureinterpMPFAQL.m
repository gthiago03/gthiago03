function [pressurinterp]=ferncodes_pressureinterpMPFAQL(p,nflagno,w,s)
global coord esurn1 esurn2


%% interpola��o das press�es nos n�s
for no=1:size(coord,1)
    nec1=esurn2(no+1)-esurn2(no);
    p1=0;
    auxflag=202; % fluxo preescrito distinto de zero, usualmente � utilizado
    % no problema Buckley-Leverett de Bastian
    if nflagno(no,1) >200
        if nflagno(no,1)==auxflag
            for j=1:nec1
                element1=esurn1(esurn2(no)+j);
                p1=p1+w(esurn2(no)+j)*p(element1);
            end
            
            % s(no,1) � contribui��o do termo fluxo distinto de zero no n�
            % "no", quando o fluxo � zero o s(no,1) � zero.
            p1=p1+s(no,1);
        else
            for j=1:nec1
                element1=esurn1(esurn2(no)+j);
                p1=p1+w(esurn2(no)+j)*p(element1);
            end
        end
        
    else
        p1=nflagno(no,2);
    end
    
    pressurinterp(no,1)=p1;
    
end

end