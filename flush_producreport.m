%--------------------------------------------------------------------------
%UNIVERSIDADE FEDERAL DE PERNAMBUCO
%CENTRO DE TECNOLOGIA E GEOCIENCIAS
%PROGRAMA DE POS GRADUACAO EM ENGENHARIA CIVIL
%TOPICOS ESPECIAIS EM DINAMICA DOS FLUIDOS COMPUTACIONAL
%--------------------------------------------------------------------------
%Subject: numerical routine to load the geometry in gmsh and to create...
%the mesh parameter 
%Type of file: FUNCTION
%Criate date: 12/08/2012
%Modify data:   /  /2012
%Advisor: Paulo Lyra and Darlan Karlo
%Programer: M�rcio Souza
%--------------------------------------------------------------------------
%Goals:
%This function writes and plots the production data.   

%--------------------------------------------------------------------------
%Additional Comments: It is called by the function "IMPES.m"


%--------------------------------------------------------------------------

function [productionreport] = flush_producreport(contime,oilflowratevec,...
    oilaccumvec,watercutvec,producelem)
%Define global parameters
global filepath resfolder;

%Initialize "productionreport" according to "producelem" size
%There is one producer well
if length(producelem) == 1
    %Initialize "productionreport"
    productionreport = zeros(length(contime) + 1,4);
    %Attribute watercut column
    productionreport(2:size(productionreport,1),4) = watercutvec;
%There are more than one producer wells
else
    productionreport = zeros(length(contime) + 1,5);
    productionreport(2:size(productionreport,1),4:5) = watercutvec(:,1:2); % quadrilateral
    % modifiquei para problema 
    %productionreport(2:size(productionreport,1),4:5) = watercutvec(:,3:4);
    % triangular
end  %End of IF

%Initialize "prodrepsize"
prodrepsize = size(productionreport,1);

%Define other parameters for "productionreport"
productionreport(2:prodrepsize,1) = contime;
productionreport(2:prodrepsize,2) = oilflowratevec;
productionreport(2:prodrepsize,3) = oilaccumvec;

%Write table (Time (VPI), Oil Flow rate, Accumulated Oil and Water Cut)
%Create the file name
prfilename = [resfolder '_' 'ProdutionReport.dat'];
%Open the file
writeproductionreport = fopen([filepath '\' prfilename],'w'); 

%Write "productionreport" according to "producelem" size
%There is one producer well
if length(producelem) == 1
    fprintf(writeproductionreport,...
        '%26.16E %26.16E %26.16E %26.16E\r\n',productionreport(:,1:4)');
%There are more than one producer wells (but just write two)
else
    fprintf(writeproductionreport,...
        '%26.16E %26.16E %26.16E %26.16E %26.16E\r\n',...
        productionreport(:,1:5)');
end  %End of IF

%Close the file "writeproductionreport.dat"
fclose(writeproductionreport);
