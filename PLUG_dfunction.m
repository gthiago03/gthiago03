
%--------------------------------------------------------------------------
%Subject: numerical routine to solve a two phase flow in porous media

%Programer: Fernando R L Contreras
%--------------------------------------------------------------------------
%Goals: this FUNCTION gets the value to boundary condition according the
%benchmark which intend to run.

%--------------------------------------------------------------------------
%Aditional comments:

%--------------------------------------------------------------------------

%Fill the matrix of permeability as a function of element to left and right
%of half edge evaluated. This function receives "kmap" and a feature of
%element which wants to know the permeability ("kfeature").
function [dmap,dmedio,gamma] = PLUG_dfunction
%Define global parameters:
global  numcase elem;

%Choose the benchmark to attribute permeability.
switch numcase
    %----------------------------------------------------------------------
    
    case 231
        
        %Build "dmap"
        dmap(1,:) = [1 12 0 0 12];
        dmedio=10;
        gamma=0;
        
    case 232
        %Build "dmap"
        dmap(1,:) = [1 0.075 0 0 0.075];
        dmedio=0.075;
        gamma=0;
    case 233
        %Build "dmap"
        
        dmap(1,:) = [1 0.5 0 0 0.5];
        dmedio=0.5;
        gamma=0;
        
    case 234
        %Build "dmap"
        
        dmap(1,:) = [1 1 0 0 1];
        dmedio=1;
        gamma=0;
    case 235
        dmap(1,:) = [1 2 0 0 2];
        dmedio=2;
        gamma=0;
    case 236
        dmap(1,:) = [1 0.05 0 0 0.05];
        dmedio=0.05;
        gamma=0.5;
    case 237
        dmap(1,:) = [1 0.05 0 0 0.05];
        dmedio=0.05;
        gamma=0.05;
    case 238
        dmap(1,:) = [1 0.05 0 0 0.05];
        dmedio=0.05;
        gamma=0.02;
    case 239
        %dmap(1,:) = [1 10 0 0 10]; % para malha menos refinada dx=0.5
        %dmedio=10;
        dmap(1,:) = [1 4 0 0 4]; % para malha mas refinada dx=0.2
        dmedio=4;
        gamma=0;
    case 241
        %dmap(1,:) = [1 4 0 0 4]; % r=1
        %dmap(1,:) = [1 4 0 0 0.4]; % r=10
        dmap(1,:) = [1 4 0 0 0.04]; % r=100
        %dmap(1,:) = [1 4 0 0 0.004]; % r=10000
        dmedio=4;
        gamma=0;
    case 242
        dmap(1,:) = [1 0.01 0 0 0.01]; % para malha menos refinada dx=0.5
        
        %dmap(1,:) = [1 0.004 0 0 0.004]; % para malha mas refinada dx=0.2
        dmedio=norm(dmap(1,2:4),inf);
        gamma=0;
    case 243
        dmap(1,:) = [1 0.01 0 0 0.01]; % para malha menos refinada dx=0.5
        dmedio=norm(dmap(1,2:4),inf);
        %dmap(1,:) = [1 0.004 0 0 0.004]; % para malha mas refinada dx=0.2
        %dmedio=0.004;
        gamma=0;
    case 244
        dmap(1,:) = [1 0.15 0 0 0.15]; % para malha menos refinada dx=0.5
        dmedio=0.15;
        %dmap(1,:) = [1 0.004 0 0 0.004]; % para malha mas refinada dx=0.2
        %dmedio=0.004;
        gamma=0;
    case 245
        for i=1:size(elem,1)
            dmap(i,:) = [i 0.01 0 0 0.01];
        end
        dmedio=norm(dmap(1,2:4),inf);
        gamma=0;
    case 247
        for i=1:size(elem,1)
            dmap(i,:) = [i 1.0000e-04 0 0 1.0000e-04];
        end
        dmedio=norm(dmap(1,2:4),inf);
        gamma=0;
    case 249
        for i=1:size(elem,1)
            dmap(i,:) = [i 1.0000e-04 0 0 1.0000e-04];
        end
        dmedio=norm(dmap(1,2:4),inf);
        gamma=0;
    case 248
        
        dmap = [1 0.01 0 0 0.01];
        dmedio=norm(dmap(1,2:4),inf);
        gamma=0;
    case 246
        for i=1:size(elem,1)
            dmap(i,:) = [i 1e-4 0 0 1e-4];
        end
        dmedio=norm(dmap(1,2:4),inf);
        gamma=0;
        
end
end
