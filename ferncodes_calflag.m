%It is called by "preMPFA.m"

function nflag = ferncodes_calflag(a)
% determinar o flag do n� interior e fronteira de Neumann
global coord bedge bcflag;
nflag = 5000*ones(size(coord,1),2);

for ifacont = 1:size(bedge,1)
    x = logical(bcflag(:,1) == bedge(ifacont,4));
    vertex = bedge(ifacont,1);
    %Second column receives the boundary condition value.
    nflag(vertex,2) = PLUG_bcfunction(vertex,x,a);
    %First column receives the boundary condition flag.
    nflag(vertex,1) = bcflag(x,1);
end  %End of FOR