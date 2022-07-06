%--------------------------------------------------------------------------
%UNIVERSIDADE FEDERAL DE PERNAMBUCO
%CENTRO DE TECNOLOGIA E GEOCIENCIAS
%PROGRAMA DE POS GRADUACAO EM ENGENHARIA CIVIL
%TOPICOS ESPECIAIS EM DINAMICA DOS FLUIDOS COMPUTACIONAL
%--------------------------------------------------------------------------
%Subject: numerical routine to solve a two phase flow in porous media
%Type of file: FUNCTION
%Criate date: 31/07/2012
%Modify data:  / /2012
%Adviser: Paulo Lyra and Darlan Karlo
%Programer: M�rcio Souza
%Modified: Fernando Contreras, 2021
%--------------------------------------------------------------------------
%Goals:

%--------------------------------------------------------------------------
%This FUNCTION calculate the

%--------------------------------------------------------------------------

function [transmvecleft,transmvecright,knownvecleft,knownvecright,storeinv,...
    Bleft,Bright,overedgecoord,normk,Fg,mapinv,maptransm,mapknownvec,...
    pointedge,bodyterm,Hesq,Kde,Kn,Kt,Ded,V,N,kmap,nflag,parameter,weightDMP,...
    nflagface,p_old,contnorm,gravelem,gravpoint,gravface] = preMPFA(kmap,klb)
%Define global parameters:
global pmethod elem coord keygravity gravresult gravrate bedge inedge flag

%Obtain the coordinate of both CENTER and AUXILARY nodes of elements which
%constitute the mash. The AREA of each element is also calculated.
%Message to user:
disp(' ');
disp('---------------------------------------------------');
disp('>> Preprocessing Pressure Equation...');
disp(' ');

%Initialize all parameters of output
transmvecleft = 0; transmvecright = 0; knownvecleft = 0; knownvecright = 0;
storeinv = 0; Bleft = 0; Bright = 0; Fg = 0; mapinv = 0; maptransm = 0;
mapknownvec = 0; pointedge = 0; bodyterm = 0; Hesq = 0; Kde = 0; Kn = 0;
Kt = 0; Ded = 0; V = 0; N = 0; nflag = 0; parameter=0; weightDMP=0;
nflagface=0;p_old=0;contnorm=0;
%Define parametric variables:
%Parameter Used in Full Pressure Support (FPS)
%"p" quadrature point to flux in the auxilary sub interaction region
p = 1;
%Parameter Used in Full Pressure Support (FPS) and Triangle Pressure
%Support (TPS)
%"q" quadrature point to flux in the sub interaction region
q = 1;

%Fill the matrix "overedgecoord"
overedgecoord = overedge;
%Define the norm of permeability tensor ("normk")
[normk,kmap] = calcnormk(kmap);

%Get the length of the edge with non-null Neumann Boundary Condition.
knownboundlength = getknownboundlength(klb);

%--------------------------------------------------------------------------
%Calculate the TRANSMISSIBILITY parameters:

%Chose the type of MPFA according "pmethod"
switch char(pmethod)
    %Calculate the transmissibilities from TPFA
    case 'tpfa'
        [Hesq,Kde,Kn,Kt,Ded] = ferncodes_Kde_Ded_Kt_Kn(kmap);
        nflag = ferncodes_calflag(0);
        %Calculate the little matrices for MPFA-TPS (Aavatsmark et al., 1998)
    case 'mpfao'
        [transmvecleft,transmvecright,knownvecleft,knownvecright,storeinv,...
            Bleft,Bright,Fg,mapinv,maptransm,mapknownvec,pointedge,...
            bodyterm] = transmTPS(kmap,q,knownboundlength);
        %Calculate the little matrices for MPFA-FPS (Edwards and Zheng, 2008)
    case 'fps'
        [transmvecleft,transmvecright,knownvecleft,knownvecright,storeinv,...
            Bleft,Bright,Fg,mapinv,maptransm,mapknownvec,pointedge,...
            bodyterm] = transmFPS(kmap,q,p,knownboundlength);
        %Calculate the little matrices for MPFA-Enriched (Chen et al., 2008)
    case 'empfa'
        [transmvecleft,transmvecright,knownvecleft,knownvecright,storeinv,...
            Bleft,Bright,Fg,mapinv,maptransm,mapknownvec,pointedge,...
            bodyterm] = transmEnriched(kmap,knownboundlength);
        %Calculate geometrical and physical terms to be used in MPFA-Diamond
        %(Gao and Wu, 2010).
    case 'mpfad'
        %Get "ferncodes_nflag"
        nflag = ferncodes_calflag;
        %Get preprocessed terms:
        [Hesq,Kde,Kn,Kt,Ded] = ferncodes_Kde_Ded_Kt_Kn(kmap);
        %Call another parameters that I don't know.
        [V,N,] = ferncodes_elementface(nflag);
    case 'mpfaql'
        % calculo dos parametros ou constantes (ksi)
        [parameter]=ferncodes_coefficient(kmap);
        % adequa��o dos flags de contorno
        nflag = ferncodes_calflag;
        % calculo dos pesos DMP
        [weightDMP]=ferncodes_weightnlfvDMP(kmap);
        %Call another parameters that I don't know.
        [V,N,] = ferncodes_elementface(nflag);
    case 'nlfvpp'
        p_old=0e1*ones(size(elem,1),1); % inicializando a press�o
        
        % utilize tolerancia menor que 10^-12 para testes de convergencia
        % Para teste de monotonicidade ou problemas bif�sicos utilize 10^-8.
                              % # itera��es de Picard
        %temos usado para muitos estes o seguinte rutina
        [parameter,contnorm]=ferncodes_coefficient(kmap);
        % adequa��o dos flags de contorno
        nflag = ferncodes_calflag(0);
        %Call another parameters that I don't know.
        [V,N,] = ferncodes_elementface(nflag);
    case 'mpfah'
        
        % faces alrededor de um elemento
        [facelement]=ferncodes_elementfacempfaH;
        %% calculoa dos pontos armonicos
        [pointarmonic]=ferncodes_harmonicopoint(kmap);
        % calculo dos parametros ou constantes (ksi)
        % temos usado este parametro durante muito tempo em muitos testes
        [parameter,auxface]=ferncodes_coefficientmpfaH(facelement,pointarmonic,kmap);
        
        % adequa��o dos flag de face de contorno
        nflagface= ferncodes_contflagface;
        % adequa��o dos nos flags de contorno
        % adequa��o dos flags de contorno
        nflag = ferncodes_calflag;
        %% calculo dos pesos DMP
        [weightDMP]=ferncodes_weightnlfvDMP(kmap);
    case 'nlfvh'
        p_old=1e1*ones(size(elem,1),1); % inicializando a press�o
        
        % utilize tolerancia menor que 10^-12 para testes de convergencia
        % Para teste de monotonicidade ou problemas bif�sicos utilize 10^-8.
                              % # itera��es de Picard
        % faces alrededor de um elemento
        [facelement]=ferncodes_elementfacempfaH;
        %% calculoa dos pontos armonicos
        [pointarmonic]=ferncodes_harmonicopoint(kmap);
        % calculo dos parametros ou constantes (ksi)
        % temos usado este parametro durante muito tempo em muitos testes
        [parameter,auxface]=ferncodes_coefficientmpfaH(facelement,pointarmonic,kmap);
        
        % adequa��o dos flag de face de contorno
        nflagface= ferncodes_contflagface;
        % adequa��o dos nos flags de contorno
        % adequa��o dos flags de contorno
        nflag = ferncodes_calflag;
        %% calculo dos pesos DMP
        [weightDMP]=ferncodes_weightnlfvDMP(kmap);
    case 'nlfvdmp'
        p_old=1e1*ones(size(elem,1),1); % inicializando a press�o
        
        % utilize tolerancia menor que 10^-12 para testes de convergencia
        % Para teste de monotonicidade ou problemas bif�sicos utilize 10^-8.
                               % # itera��es de Picard
        % faces alrededor de um elemento
        [facelement]=ferncodes_elementfacempfaH;
        %% calculoa dos pontos armonicos
        [pointarmonic]=ferncodes_harmonicopoint(kmap);
        % calculo dos parametros ou constantes (ksi)
        % temos usado este parametro durante muito tempo em muitos testes
        [parameter,auxface]=ferncodes_coefficientmpfaH(facelement,pointarmonic,kmap);
        
        % adequa��o dos flag de face de contorno
        nflagface= ferncodes_contflagface;
        % adequa��o dos nos flags de contorno
        % adequa��o dos flags de contorno
        nflag = ferncodes_calflag(0);
        %% calculo dos pesos DMP
        [weightDMP]=ferncodes_weightnlfvDMP(kmap);
        
        if strcmp(keygravity,'n')
            gravresult=zeros(size(elem,1),1);
            gravrate=zeros(size(bedge,1)+size(inedge,1),1);
            gravelem=zeros(size(elem,1),1);
            gravpoint=zeros(size(coord,1),1);
            gravface=zeros(size(bedge,1)+size(inedge,1),1);
            
        else
            % Thiago
            %Gt = jtgs_gravitationalterm(overedgecoord, kmap);
            % Fernando
            [vec_gravelem,vec_gravface,vec_gravpoint,gravelem,gravpoint,gravface]=PLUG_Gfunction;
            [gravresult, gravrate]=gravitation(kmap,vec_gravelem,vec_gravface);
            [nflagface] = standard_discr(kmap,vec_gravelem,vec_gravface,vec_gravpoint,nflagface);
            flag = nflagface;
            %Gt=-Gt;
            nflag(5,2) = (nflagface(1,2)+nflagface(2,2))/2;
            nflag(6,2) = (nflagface(2,2)+nflagface(3,2))/2;
            nflag(7,2) = (nflagface(3,2)+nflagface(4,2))/2;
            nflag(11,2) = (nflagface(9,2)+nflagface(10,2))/2;
            nflag(12,2) = (nflagface(10,2)+nflagface(11,2))/2;
            nflag(13,2) = (nflagface(11,2)+nflagface(12,2))/2;
            
        end
        
end  %End of SWITCH

%Message to user:
disp('>> "preMPFA" was finished with success!');

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%FUNCTION DEFINITION
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%FUNCTION "overedgesection"
%--------------------------------------------------------------------------

%This function find the half point of each straight line. After that obtain
%%its coordinate.
%"nodeval" is the node evaluated. This is one of two components calculated
%by "interselemnode"
function [coordsection] = overedgesection(edgematrix)
%Define global parameters:
global coord;

%Initialize the matrix used in this function
coordsection = zeros(size(edgematrix,1),3);
%This loop swept among the limits stabelished above
%"firstlim" and "lastlim" are parameters which define where the loop begins
%and where its finish

iover = 1:size(edgematrix,1);
coordsection(iover,1:3) = 0.5*(coord(edgematrix(iover,1),:) + ...
    coord(edgematrix(iover,2),:));

%--------------------------------------------------------------------------
%Function "overedge"
%--------------------------------------------------------------------------

function [overedgecoord] = overedge
%Define global parameters:
global bedge inedge;

%Initialize the matrix "overedgecoord"
overedgecoord = zeros((size(bedge,1) + size(inedge,1)),3);

%Fill "overedgecoord" (just edge over boundary)
overedgecoord(1:size(bedge,1),:) = overedgesection(bedge);
%Fill the "overedgecoord" rest (just edge inside domain)
%"continedge" is an internal edge's counter
overedgecoord(size(bedge,1) + 1:size(overedgecoord,1),:) = ...
    overedgesection(inedge);

%--------------------------------------------------------------------------
%Function "calcnormk"
%--------------------------------------------------------------------------

function [normk,kmap] = calcnormk(kmap)
%Define global parameters:
global elem centelem;

%Initialize "normk" (it is a vector)
normk = zeros(size(centelem,1),1);
%Define the norm of permeability tensor
%Obtain "kmap" for each case
kmap = PLUG_kfunction(kmap);
%Swept all elements
for ik = 1:length(normk)
    %Define the material pointer in "elem"
    pointer = elem(ik,5);
    %It catches only the permeability components
    permcompon = [kmap(pointer,2) kmap(pointer,3); ...
        kmap(pointer,4) kmap(pointer,5)];
    %Calculate the norm of tensor
    normk(ik) = norm(permcompon);
end  %End of FOR
