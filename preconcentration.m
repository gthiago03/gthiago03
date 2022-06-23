
%Modified: Fernando Contreras, 2021


function [bedge,bcflagc,wells]= preconcentration(bedge,wells)
global  coord numcase


%% flags adicionais para o problema 5.3.1-5.3.2
switch numcase
    case 231
        bcflagc(1,1)=251;
        bcflagc(1,2)=0;
        bcflagc(2,1)=50;  % dicichlet concen
        bcflagc(2,2)=10; % valor dirich
        bcflagc(3,1)=250; % dirich neumann
        bcflagc(3,2)=0;  % valor concent
        
        c=logical(bedge(:,4)==201);
        bedge([c],6)=251;
        d=logical(bedge(:,4)==101);
        bedge([d],6)=50;
        e=logical(bedge(:,4)==102);
        bedge([e],6)=250;
        
        c1=logical(bedge(:,5)==201);
        bedge([c1],7)=251;
        d1=logical(bedge(:,5)==101);
        bedge([d1],7)=50;
        e1=logical(bedge(:,5)==102);
        bedge([e1],7)=250;
        wells=0;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % bcflag(6,1)=50;  % dicichlet concen
        % bcflag(7,1)=51; % dirich concen
        % bcflag(8,1)=250;  % dicichlet concen
        % bcflag(9,1)=251; % dirich concen
        % bcflag(6,2)=10; % valor dirich
        % bcflag(7,2)=10;  % valor concent
        % bcflag(8,2)=0; % valor dirich
        % bcflag(9,2)=0;  % valor concent
        %
        % c=logical(bedge(:,4)==201);
        % bedge([c],6)=201;
        % d=logical(bedge(:,4)==101);
        % bedge([d],6)=50;
        % e=logical(bedge(:,4)==102);
        % bedge([e],6)=51;
        % f=logical(bedge(:,4)==103);
        % bedge([f],6)=250;
        % g=logical(bedge(:,4)==104);
        % bedge([g],6)=251;
        %
        % c1=logical(bedge(:,5)==201);
        % bedge([c1],7)=201;
        % d1=logical(bedge(:,5)==101);
        % bedge([d1],7)=50;
        % e1=logical(bedge(:,5)==102);
        % bedge([e1],7)=51;
        % f1=logical(bedge(:,5)==103);
        % bedge([f1],7)=250;
        % g1=logical(bedge(:,5)==104);
        % bedge([g1],7)=251;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 232
        bcflagc(1,1)=251;
        bcflagc(1,2)=0;
        bcflagc(2,1)=50;  % dicichlet concen
        bcflagc(2,2)=10; % valor dirich
        bcflagc(3,1)=250; % dirich neumann
        bcflagc(3,2)=0;  % valor concent
        
        c=logical(bedge(:,4)==201);
        bedge([c],6)=251;
        d=logical(bedge(:,4)==101);
        bedge([d],6)=50;
        e=logical(bedge(:,4)==102);
        bedge([e],6)=250;
        
        c1=logical(bedge(:,5)==201);
        bedge([c1],7)=251;
        d1=logical(bedge(:,5)==101);
        bedge([d1],7)=50;
        e1=logical(bedge(:,5)==102);
        bedge([e1],7)=250;
        wells=0;
    case 233
        bcflagc(1,1)=251;
        bcflagc(1,2)=0;
        bcflagc(2,1)=50;  % dicichlet concen
        bcflagc(2,2)=10; % valor dirich
        bcflagc(3,1)=250; % dirich neumann
        bcflagc(3,2)=0;  % valor concent
        
        c=logical(bedge(:,4)==201);
        bedge([c],6)=251;
        d=logical(bedge(:,4)==101);
        bedge([d],6)=50;
        e=logical(bedge(:,4)==102);
        bedge([e],6)=250;
        
        c1=logical(bedge(:,5)==201);
        bedge([c1],7)=251;
        d1=logical(bedge(:,5)==101);
        bedge([d1],7)=50;
        e1=logical(bedge(:,5)==102);
        bedge([e1],7)=250;
        wells=0;
        
    case 234
        bcflagc(1,1)=251;
        bcflagc(1,2)=0;
        bcflagc(2,1)=50;  % dicichlet concen
        bcflagc(2,2)=10; % valor dirich
        bcflagc(3,1)=250; % dirich neumann
        bcflagc(3,2)=0;  % valor concent
        
        c=logical(bedge(:,4)==201);
        bedge([c],6)=251;
        d=logical(bedge(:,4)==101);
        bedge([d],6)=50;
        e=logical(bedge(:,4)==102);
        bedge([e],6)=250;
        
        c1=logical(bedge(:,5)==201);
        bedge([c1],7)=251;
        d1=logical(bedge(:,5)==101);
        bedge([d1],7)=50;
        e1=logical(bedge(:,5)==102);
        bedge([e1],7)=250;
        wells=0;
    case 235
        bcflagc(1,1)=251;
        bcflagc(1,2)=0;
        bcflagc(2,1)=50;  % dicichlet concen
        bcflagc(2,2)=10; % valor dirich
        bcflagc(3,1)=250; % dirich neumann
        bcflagc(3,2)=0;  % valor concent
        
        c=logical(bedge(:,4)==201);
        bedge([c],6)=251;
        d=logical(bedge(:,4)==101);
        bedge([d],6)=50;
        e=logical(bedge(:,4)==102);
        bedge([e],6)=250;
        
        c1=logical(bedge(:,5)==201);
        bedge([c1],7)=251;
        d1=logical(bedge(:,5)==101);
        bedge([d1],7)=50;
        e1=logical(bedge(:,5)==102);
        bedge([e1],7)=250;
        wells=0;
    case 236
        bcflagc(1,1)=251;
        bcflagc(1,2)=0;
        bcflagc(2,1)=50;  % dicichlet concen
        bcflagc(2,2)=10; % valor dirich
        bcflagc(3,1)=250; % dirich neumann
        bcflagc(3,2)=0;  % valor concent
        
        c=logical(bedge(:,4)==201);
        bedge([c],6)=251;
        d=logical(bedge(:,4)==101);
        bedge([d],6)=50;
        e=logical(bedge(:,4)==102);
        bedge([e],6)=250;
        
        c1=logical(bedge(:,5)==201);
        bedge([c1],7)=251;
        d1=logical(bedge(:,5)==101);
        bedge([d1],7)=50;
        e1=logical(bedge(:,5)==102);
        bedge([e1],7)=250;
        wells=0;
    case 237
        bcflagc(1,1)=251;
        bcflagc(1,2)=0;
        bcflagc(2,1)=50;  % dicichlet concen
        bcflagc(2,2)=10; % valor dirich
        bcflagc(3,1)=250; % dirich neumann
        bcflagc(3,2)=0;  % valor concent
        
        c=logical(bedge(:,4)==201);
        bedge([c],6)=251;
        d=logical(bedge(:,4)==101);
        bedge([d],6)=50;
        e=logical(bedge(:,4)==102);
        bedge([e],6)=250;
        
        c1=logical(bedge(:,5)==201);
        bedge([c1],7)=251;
        d1=logical(bedge(:,5)==101);
        bedge([d1],7)=50;
        e1=logical(bedge(:,5)==102);
        bedge([e1],7)=250;
        wells=0;
    case 238
        bcflagc(1,1)=251;
        bcflagc(1,2)=0;
        bcflagc(2,1)=50;  % dicichlet concen
        bcflagc(2,2)=10; % valor dirich
        bcflagc(3,1)=250; % dirich neumann
        bcflagc(3,2)=0;  % valor concent
        
        c=logical(bedge(:,4)==201);
        bedge([c],6)=251;
        d=logical(bedge(:,4)==101);
        bedge([d],6)=50;
        e=logical(bedge(:,4)==102);
        bedge([e],6)=250;
        
        c1=logical(bedge(:,5)==201);
        bedge([c1],7)=251;
        d1=logical(bedge(:,5)==101);
        bedge([d1],7)=50;
        e1=logical(bedge(:,5)==102);
        bedge([e1],7)=250;
        wells=0;
    case 239
        bcflagc(1,1)=50;% dicichlet concen
        bcflagc(1,2)=0;% valor dirich
        bcflagc(2,1)=51; % dicichlet concen 
        bcflagc(2,2)=1; % valor dirich
        bcflagc(3,1)=250; % neumann boundary
        bcflagc(3,2)=0;  % valor 
        bcflagc(4,1)=251;  % neumann boundary
        bcflagc(4,2)=0; % valor 
        
        c=logical(bedge(:,4)==101);
        bedge([c],6)=50;
        d=logical(bedge(:,4)==102);
        bedge([d],6)=51;
        e=logical(bedge(:,4)==103);
        bedge([e],6)=250;
        f=logical(bedge(:,4)==201);
        bedge([f],6)=251;
        
        c1=logical(bedge(:,5)==101);
        bedge([c1],7)=50;
        d1=logical(bedge(:,5)==102);
        bedge([d1],7)=51;
        e1=logical(bedge(:,5)==103);
        bedge([e1],7)=250;
        f1=logical(bedge(:,5)==201);
        bedge([f1],7)=251;
        wells=0;
    case 241
        bcflagc(1,1)=50;% dicichlet concen
        bcflagc(1,2)=0;% valor dirich
        bcflagc(2,1)=51; % dicichlet concen 
        bcflagc(2,2)=1; % valor dirich
        bcflagc(3,1)=250; % neumann boundary
        bcflagc(3,2)=0;  % valor 
        bcflagc(4,1)=251;  % neumann boundary
        bcflagc(4,2)=0; % valor 
        
        c=logical(bedge(:,4)==101);
        bedge([c],6)=50;
        d=logical(bedge(:,4)==102);
        bedge([d],6)=51;
        e=logical(bedge(:,4)==103);
        bedge([e],6)=250;
        f=logical(bedge(:,4)==201);
        bedge([f],6)=251;
        
        c1=logical(bedge(:,5)==101);
        bedge([c1],7)=50;
        d1=logical(bedge(:,5)==102);
        bedge([d1],7)=51;
        e1=logical(bedge(:,5)==103);
        bedge([e1],7)=250;
        f1=logical(bedge(:,5)==201);
        bedge([f1],7)=251;
        wells=0;
    case 242
        bcflagc(1,1)=50;% dicichlet concen
        bcflagc(1,2)=0;% valor dirich
        bcflagc(2,1)=51; % dicichlet concen 
        bcflagc(2,2)=1; % valor dirich
        bcflagc(3,1)=250; % neumann boundary
        bcflagc(3,2)=0;  % valor 
        bcflagc(4,1)=251;  % neumann boundary
        bcflagc(4,2)=0; % valor 
        
        c=logical(bedge(:,4)==101);
        bedge([c],6)=50;
        d=logical(bedge(:,4)==102);
        bedge([d],6)=51;
        e=logical(bedge(:,4)==103);
        bedge([e],6)=250;
        f=logical(bedge(:,4)==201);
        bedge([f],6)=251;
        
        c1=logical(bedge(:,5)==101);
        bedge([c1],7)=50;
        d1=logical(bedge(:,5)==102);
        bedge([d1],7)=51;
        e1=logical(bedge(:,5)==103);
        bedge([e1],7)=250;
        f1=logical(bedge(:,5)==201);
        bedge([f1],7)=251;
        wells=0;
    case 243
        bcflagc(1,1)=50;% dicichlet concen
        bcflagc(1,2)=0;% valor dirich
        bcflagc(2,1)=51; % dicichlet concen 
        bcflagc(2,2)=1; % valor dirich
        bcflagc(3,1)=250; % neumann boundary
        bcflagc(3,2)=0;  % valor 
        bcflagc(4,1)=251;  % neumann boundary
        bcflagc(4,2)=0; % valor 
        
        c=logical(bedge(:,4)==101);
        bedge([c],6)=50;
        d=logical(bedge(:,4)==102);
        bedge([d],6)=51;
        e=logical(bedge(:,4)==103);
        bedge([e],6)=250;
        f=logical(bedge(:,4)==201);
        bedge([f],6)=251;
        
        c1=logical(bedge(:,5)==101);
        bedge([c1],7)=50;
        d1=logical(bedge(:,5)==102);
        bedge([d1],7)=51;
        e1=logical(bedge(:,5)==103);
        bedge([e1],7)=250;
        f1=logical(bedge(:,5)==201);
        bedge([f1],7)=251;
%         bcflagc(1,1)=251;
%         bcflagc(1,2)=0;
%         bcflagc(2,1)=50;  % dicichlet concen
%         bcflagc(2,2)=10; % valor dirich
%         bcflagc(3,1)=250; % dirich neumann
%         bcflagc(3,2)=0;  % valor concent
%         
%         c=logical(bedge(:,4)==201);
%         bedge([c],6)=251;
%         d=logical(bedge(:,4)==101);
%         bedge([d],6)=50;
%         e=logical(bedge(:,4)==102);
%         bedge([e],6)=250;
%         
%         c1=logical(bedge(:,5)==201);
%         bedge([c1],7)=251;
%         d1=logical(bedge(:,5)==101);
%         bedge([d1],7)=50;
%         e1=logical(bedge(:,5)==102);
%         bedge([e1],7)=250;
        wells=0;
    case 244
        bcflagc(1,1)=250;
        bcflagc(1,2)=0;
        
        c=logical(bedge(:,4)==101);
        bedge([c],6)=250;
        
        c1=logical(bedge(:,5)==101);
        bedge([c1],7)=250;
        wells (1,1:6)=[1201 1 301 1 0 1];
    case 245
        bcflagc(1,1)=250;
        bcflagc(1,2)=0;
         c=logical(bedge(:,4)==201);
        bedge([c],6)=250;
        c1=logical(bedge(:,5)==201);
        bedge([c1],7)=250;
   case 247
        bcflagc(1,1)=250;
        bcflagc(1,2)=0;
         c=logical(bedge(:,4)==201);
        bedge([c],6)=250;
        c1=logical(bedge(:,5)==201);
        bedge([c1],7)=250;
   case 248
        bcflagc(1,1)=50;
        bcflagc(1,2)=0;
         c=logical(bedge(:,4)==101);
        bedge([c],6)=50;
        c1=logical(bedge(:,5)==101);
        bedge([c1],7)=50;
    case 249
        
         
        bcflagc(1,1)=251;
        bcflagc(1,2)=0;
        bcflagc(2,1)=50;  % dicichlet concen
        bcflagc(2,2)=1; % valor dirich
        bcflagc(3,1)=250; % dirich neumann
        bcflagc(3,2)=0;  % valor concent
        
        c=logical(bedge(:,4)==201);
        bedge([c],6)=251;
        d=logical(bedge(:,4)==101);
        bedge([d],6)=50;
        e=logical(bedge(:,4)==102);
        bedge([e],6)=250;
        
        c1=logical(bedge(:,5)==201);
        bedge([c1],7)=251;
        d1=logical(bedge(:,5)==101);
        bedge([d1],7)=50;
        e1=logical(bedge(:,5)==102);
        bedge([e1],7)=250;
        j=1;
        for i=1:size(bedge,1)
           if bedge(i,7)==250 
              wells(j,1)=bedge(i,3);
              wells(j,2)=2;
              wells(j,3)=0;
              wells(j,4)=0;
              wells(j,5)=501;
              wells(j,6)=0;
              j=j+1;
           end            
        end
   case 246
        bcflagc(1,1)=50;% dicichlet concen
        bcflagc(1,2)=0;% valor dirich
        bcflagc(2,1)=51; % dicichlet concen 
        bcflagc(2,2)=1; % valor dirich
        bcflagc(3,1)=250; % neumann boundary
        bcflagc(3,2)=0;  % valor 
        bcflagc(4,1)=251;  % neumann boundary
        bcflagc(4,2)=0; % valor 
        
        c=logical(bedge(:,4)==101);
        bedge([c],6)=50;
        d=logical(bedge(:,4)==102);
        bedge([d],6)=51;
        e=logical(bedge(:,4)==103);
        bedge([e],6)=250;
        f=logical(bedge(:,4)==201);
        bedge([f],6)=251;
        
        c1=logical(bedge(:,5)==101);
        bedge([c1],7)=50;
        d1=logical(bedge(:,5)==102);
        bedge([d1],7)=51;
        e1=logical(bedge(:,5)==103);
        bedge([e1],7)=250;
        f1=logical(bedge(:,5)==201);
        bedge([f1],7)=251;
end
% end flags
%% flags adicionais para o problema 5.3.3 Darlan
% if strcmp(Prob,'2D')
%     bcflag(8,1)=50; % flag aux para a concentrac�o no m�dio
%     bcflag(9,1)=51; % flag aux para a concentrac�o acima e embaixo
%     bcflag(10,1)=52; % flag aux para a concentrac�o acima e embaixo
%     bcflag(11,1)=250; % flag aux para a concentrac�o face 2
%     bcflag(12,1)=251; % flag aux para a concentrac�o face 2
%     bcflag(13,1)=252; % flag aux para a concentrac�o face 2
%     bcflag(8,2)=0; % valor da Concentra��o em bcflag(9,1)
%     bcflag(9,2)=0; % valor da Concentra��o em bcflag(10,1)
%
%     % if t_old < 0.2
%     bcflag(10,2)=1; % valor da Concentra��o em bcflag(11,1)
%     %  else
%     %      bcflag(10,2)=0; % valor da Concentra��o em bcflag(11,1)
%     %  end
%     bcflag(11,2)=0; % valor da Concentra��o em bcflag(11,1)
%     bcflag(12,2)=0; % valor da Concentra��o em bcflag(11,1)
%     bcflag(13,2)=0; % valor da Concentra��o em bcflag(11,1)
%
%     c=logical(bedge(:,4)==101);
%     bedge([c],6)=50;
%     d=logical(bedge(:,4)==102);
%     bedge([d],6)=51;
%     e=logical(bedge(:,4)==103);
%     bedge([e],6)=52;
%     f=logical(bedge(:,4)==104);
%     bedge([f],6)=250;
%     g=logical(bedge(:,4)==105);
%     bedge([g],6)=251;
%     h=logical(bedge(:,4)==106);
%     bedge([h],6)=252;
%     i=logical(bedge(:,4)==201);
%     bedge([i],6)=201;
%
%     c1=logical(bedge(:,5)==101);
%     bedge([c1],7)=50;
%     d1=logical(bedge(:,5)==102);
%     bedge([d1],7)=51;
%     e1=logical(bedge(:,5)==103);
%     bedge([e1],7)=52;
%     f1=logical(bedge(:,5)==104);
%     bedge([f1],7)=250;
%     g1=logical(bedge(:,5)==105);
%     bedge([g1],7)=251;
%     h1=logical(bedge(:,5)==106);
%     bedge([h1],7)=252;
%     i1=logical(bedge(:,5)==201);
%     bedge([i1],7)=201;
% end
% % end flags
% %%
%
% nflagp=5000*ones(size(coord,1),2); % flags 7.15934 � uma constante
% nflagc=5000*ones(size(coord,1),2); % flags 7.15934 � uma constante
%
% for ifacont=1:size(bedge,1)
%
%     x=bcflag(:,1)==bedge(ifacont,4);
%     r=find(x==1);
%     nflagp(bedge(ifacont,1),2)=bcflag(r,2);
%     nflagp(bedge(ifacont,1),1)=bcflag(r,1);
%
%     x=bcflag(:,1)==bedge(ifacont,6);
%     r=find(x==1);
%     nflagc(bedge(ifacont,1),2)=bcflag(r,2);
%     nflagc(bedge(ifacont,1),1)=bcflag(r,1);
% end
end