%--------------------------------------------------------------------------
%UNIVERSIDADE FEDERAL DE PERNAMBUCO
%CENTRO DE TECNOLOGIA E GEOCIENCIAS
%PROGRAMA DE POS GRADUACAO EM ENGENHARIA CIVIL
%TOPICOS ESPECIAIS EM DINAMICA DOS FLUIDOS COMPUTACIONAL
%--------------------------------------------------------------------------
%Subject: numerical code used to simulate fluid flow in porous media. That 
%routine calls several others which defines how the equation will be solved 
%Type of file: MAIN
%Criate date: 29/02/2012
%Modify data:   /  /2012
%Advisor: Paulo Lyra and Darlan Karlo
%Programer: Márcio Souza
%Modified: Fernando Contreras, 2021
%--------------------------------------------------------------------------
%Goals: Do the manegement of simulator. This is a MAIN program. 

%--------------------------------------------------------------------------
%In this numerical routine the flow may be simulated with one or two phase 
%(water-oil generaly). The functions below are organized in order give 
%flexibility to software resourcers.
%For example: the saturation and pressure fields are calculated by IMPES,
%but it may be calculated also by a fully implicit scheme. This change is
%done in the present rountine just call each function.
%--------------------------------------------------------------------------

%Clear the screem 
clc;
%Clear all the memory of the matlab
clear all;
%Define the format of out data
format short;
%It begins the time counter and "profile".
tic
% profile on -history;

%--------------------------------------------------------------------------
%Define the global variables:
global coord centelem elem centface esurn1 esurn2 nsurn1 nsurn2 bedge inedge ...
    normals un_normals esureface1 esureface2 esurefull1 esurefull2 elemarea dens ...
    visc satlimit pormap bcflag courant totaltime filepath resfolder ...
    numcase pmethod smethod phasekey order timeorder auxcvfactor ...
    interptype multdopt goefreeopt lsneightype lsexp ...
    recovtype keygravity g benchkey rowposit nltol maxiter acel bcflagc...
    P_old P_new cpres Gt A_old ep eq p_exact flowrate flowresult u0 gravrate gravresult;

%--------------------------------------------------------------------------
%Call the "preprocessor" function

%This function obtains from gmsh file all data structure.
[coord,centelem,centface,elem,esurn1,esurn2,nsurn1,nsurn2,bedge,inedge,normals, un_normals...
    esureface1,esureface2,esurefull1,esurefull2,elemarea,dens,visc,...
    satlimit,pormap,bcflag,courant,totaltime,numcase,phasekey,pmethod,...
    smethod,xyrz,r0,symaxe,keymsfv,coarseratio,auxcvfactor,interptype,...
   multdopt,goefreeopt,order,timeorder,recovtype,lsneightype,...
    lsexp,keygravity,g,keycapil,ncaplcorey,filepath,resfolder,benchkey,...
    kmap,wells,klb,limiterflag,rowposit,nltol,maxiter,acel] = preprocessor;
if numcase==23.3 
% x=bedge(71:78,1); % UTILIZE Benchmark23_3_18_18.msh
% y=bedge(71:78,2);
% bedge(71:78,1)=y;
% bedge(71:78,2)=x;
% bedge(71:78,4:5)=102; % 18x18
% bcflag(2,1)=102;
% bcflag(2,2)=2;

%=====================================

%  x=bedge(135:150,1); % UTILIZE Benchmark23_3_FINA36.msh
%  y=bedge(135:150,2);
%  bedge(135:150,1)=y;
%  bedge(135:150,2)=x;
%  bedge(135:150,4:5)=102; % 36x36
%  bcflag(2,1)=102;
%  bcflag(2,2)=2;
%===============================================================

% unstructured mesh com furo reordenando o sentido da fronterira no
% contorno interior
 x=bedge(289:320,1);
 y=bedge(289:320,2);
 bedge(289:320,1)=y;
 bedge(289:320,2)=x;
 bedge(289:320,4:5)=102; % benchmark23_3 72x72
 bcflag(2,1)=102;
 bcflag(2,2)=2;
end
%adeSPE; % para um campo de permeabilidade da SPE active descomente.
%--------------------------------------------------------------------------
%Call "setmethod"
%elem(:,5)=1;
%Initialize general parameters:
elemsize = size(elem,1);
bedgesize = size(bedge,1);
inedgesize = size(inedge,1);

%p_exact = jtgs_analytical();
u0=1; %pressão de referência
%It preprocess the schemes and set a One-phase or Two-phase simulation.
setmethod(kmap,wells,'i',8,limiterflag,klb,elemsize,bedgesize,inedgesize);

