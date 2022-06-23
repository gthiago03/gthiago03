%--------------------------------------------------------------------------
%UNIVERSIDADE FEDERAL DE PERNAMBUCO
%CENTRO DE TECNOLOGIA E GEOCIENCIAS
%PROGRAMA DE POS GRADUACAO EM ENGENHARIA CIVIL
%TOPICOS ESPECIAIS EM DINAMICA DOS FLUIDOS COMPUTACIONAL
%--------------------------------------------------------------------------
%Subject: numerical routine to load the geometry in gmsh and to create...
%the mesh parameter 
%Type of file: FUNCTION
%Criate date: 04/12/2013
%Modify data:   /  /2013
%Advisor: Paulo Lyra and Darlan Karlo
%Programer: M�rcio Souza
%Modified: Fernando Contreras, 2021
%--------------------------------------------------------------------------
%Goals:
%    

%--------------------------------------------------------------------------
%Additional comments: the "intertype" options are:
%0. Return only the known values
%1. Volume Average;
%2. Shared Volume Average;
%3. Linear Interpolation.

%--------------------------------------------------------------------------

function [satonvertices,satonedges,flagknownvert,flagknownedge] = ...
    getsatandflag(satinbound,injecelem,Sw,nflagc,nflagfacec)
%Define global parameters:
global coord bcflagc bedge inedge elemarea; 

%Initialize "bedgesize" and "inedgesize"
bedgesize = size(bedge,1);
inedgesize = size(inedge,1);
coordsize = size(coord,1);

%Initialize the vectors "satonvertices" and "satonedges"
satonvertices = nflagc(:,2);

satonedges = zeros(bedgesize + inedgesize,1);
%Initialize "knownvalinvert" and "knownvalinedge". It indicate if there is 
%or not a prescribed value.
knownvalinvert =  nflagc(:,3);
knownvalinedge = zeros(bedgesize,1);

%Evaluate if there is non-null Neumann boundary condition (non-null Neumann 
%flux, in pressure indicate saturation, by Dirichlet boundary condition 
%prescribed).
pointneumann = logical(bcflagc(:,1) > 200 & bcflagc(:,1) < 300 & ...
        bcflagc(:,2) ~= 0);

%--------------------------
%Saturation on the VERTICES

%Verify if some vertex or edge is submited to Dirichlet boundary condition
%(in SATURATION equation).
if any(pointneumann) && any(bcflagc(pointneumann,2) > 0) && ...
        (any(nflagfacec) || length(nflagfacec) > 1)
    %Find the flags for non-null Neumann boundary condition.
    nonullflag = bcflag(pointneumann,1);

    %Swept all the edges on the boundary
    for i = 1:bedgesize
        %The element belongs to vector "injecelem" and the flag of the 5th 
        %column of "bedge" matches to "nonullflag"
        if (ismember(bedge(i,3),injecelem)) && (bedge(i,5) == nonullflag)
            %It points the position in "injecelem" corresponding the
            %element evaluated.
            pointposit = logical(injecelem == bedge(i,3));
            %Attribute to "satonedge" the known value of saturation.
            satonedges(i) = satinbound(pointposit);
            
            %It is an indication there is a known value over edge
            knownvalinedge(i) = 1;
            
            %Attribute the boundary condition for the vertices
            satonvertices(bedge(i,1:2)) = satinbound(pointposit);
            %It is an indication there is a known value over vertex
            knownvalinvert(bedge(i,1:2)) = 1;
        end  %End of IF            
    end  %End of FOR
end  %End of IF

%Get the edges whose boudary condition were NOT attributed.
iedge = 1:bedgesize;
unknownedge = iedge(logical(knownvalinedge == 0));    
%Get the vertices whose boudary condition were NOT attributed.
ivert = 1:coordsize;
unknownvert = ivert(logical(knownvalinvert == 0));    


    %Calculate the saturation into each unknown vertex    
    for i = 1:length(unknownvert)
        inode = unknownvert(i);
        %get the amount of elements surrounding the vertex evaluated.
        [esurn,] = getsurnode(inode); 
        %get the saturation into vertex evaluated.
        satinvertex = getsatinvertex(elemarea(esurn),Sw(esurn));
        %Attribute the value calculated to "satonvertices"
        satonvertices(inode) = satinvertex;
    end  %End of FOR (swept the vertices)    


%-----------------------
%Saturation on the EDGES

%Calculate the saturation into each unknown edge (just in "bedge")    
%The Saturation on the Midedge is obtained by Shared Volume Average.

    %Swept "bedge"
    for i = 1:length(unknownedge)
        %Calculate the saturation into each unknown edge ("bedge")    
        %satonedges(unknownedge(i)) = nflagfacec(i,2);
         satonedges(unknownedge(i))= 0.5*(satonvertices(bedge(i,1)) + satonvertices(bedge(i,2)));
        %It is an indication there is a known value over edge
        if nflagfacec(i,1)<100
        knownvalinedge(i)=1;
        else
        knownvalinedge(i)=0;    
        end
    end  %End of FOR (swept the edges from "bedge")

    %Swept "inedge"
    for i = 1:inedgesize
        %Calculate the saturation into each unknown edge ("inedge")    
        satonedges(bedgesize + i) = ...
            0.5*(satonvertices(inedge(i,1)) + satonvertices(inedge(i,2)));
        
    end  %End of FOR (swept the edges from "inedge")    


%Fill the flag vectors:
flagknownvert = knownvalinvert;
flagknownedge = knownvalinedge;

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%FUNCTION DEFINITION:
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%Function "getsativertex"
%--------------------------------------------------------------------------

function [satinvertex] = getsatinvertex(esurnarea,Sw_esurn)
%It uses a mean weighted by volumes
satinvertex = Sw_esurn'*esurnarea/sum(esurnarea);
