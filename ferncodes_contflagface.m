function [nflag]= ferncodes_contflagface
 
% determinar o flag do n� interior e fronteira de Neumann
global coord bedge bcflag;
nflag = 5000*ones(size(coord,1),2);

for ifacont = 1:size(bedge,1)
    x = logical(bcflag(:,1) == bedge(ifacont,5));
    vertex = bedge(ifacont,1:2);
    %Second column receives the boundary condition value.
    nflag(ifacont,2) = PLUG_bcfunction(vertex,x);
    %First column receives the boundary condition flag.
    nflag(ifacont,1) = bcflag(x,1);
end  %End of FOR
end