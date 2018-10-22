%In this script, we will analyze the use of preconditioners 
%on our 2D and 3D heat conducting system.  

%% Paths and clearing 
clear all %If we want to reset

mydir  = pwd;
idcs   = strfind(mydir,'/');
newdir = mydir(1:idcs(end)-1); %navigate to folder above
addpath([newdir '/codeTask1']); %add appropriate directory 
addpath Utils
%for previous task
%% GENERATING A NODAL MATRIX
dim_flag = 2;
%Dimension of the brain that we want to creat a nodal matrix fore
%1 => 1D, $2 => 2D, $3 => 3D

visualize = 0; %If we want to visualize or results 

switch(dim_flag)
    case 1
        %Generate a 1D nodal Matrix
    case 2
        %Generate a 2D nodal Matrix
        Nx= 100; %Nodes in x direction
        Ny= 100; %Nodes in y direction
        dx = 1; %Dicretization in x direction
        dy = 1; %Dicretization in y direction
        
        % 2D Heat Process
        tmsize = 15; %Dimension of tumor square in center
        
        [p,x_start,t_start,t_stop,max_dt_FE] = ...
            getParam_2DHeating(dx,dy,Nx,Ny,tmsize);
        
        A = p.A; %our nodal matrix
        b = p.B; %single excitation in the center
    case 3
        %Generate a 3D nodal Matrix
        N = 22; %Dimension in all 3, 3D dimensions
        R = zeros(3,N);
        R(1,:) = 1;
        R(2,:) = 0:N-1;
        R(3,:) = 1:N; 
        %assume all resistors are connected by value 1 in a line
        
        [E,~,alpha] = constructEIs(R,[],N); 
        G = E*alpha*E';
        %construct the 1D nodal Matrix 
        
        I = eye(N);
        G2 = kron(I,G) + kron(G,I);
        %construct the 2D nodal Matrix
        
        I2 = eye(N^2);
        A = kron(I2,G) + kron(G2,I);
        %construct the 3D Nodal Matrix 
        %for now, generate generic 3D nodal matrix just to test
        %preconditioner 
        
        b = zeros(N^3,1);
        b(floor(N^2/2),1) = 1; 
        %apply heat to the middle of the outermost square
        %like applying heat to the edge of the brain
        
end

%% TGCR BASELINE SOLVE
tol = 1e-6;
maxiters = 1000;

tic
[x,r_norms] = tgcr(A,b,tol,maxiters);
tgcr_time = toc;

if(visualize)
    figure
    switch(dim_flag)
        case 1
            %visualize the result in 1d
        case 2    
            imagesc(reshape(x,Nx,Ny)) 
            %visualize result of 2D
        case 3
            plot(x) 
            %visuzlize the result in 3D
    end
end
%% TGCR Solution with Preconditioner
for ii = 1:3 %Test each preconditioner
    precon_flag = ii;
    %indicates which preconditioner we wish to use
    %1 => Jacobi preconditioner 
    %2 => tridiag preconditioner 
    %3 => lowerdiag preconditioner 

    A_hat = genPrecon(precon_flag,A);

    %solve system with preconditioner
    tic
    [x_p,r_norms_p] = tgcr_P(A,b,tol,maxiters,A_hat);
    precon_time(ii) = toc;
    
    if(visualize)
        figure
        switch(dim_flag)
            case 1
                %visualize the result in 1d
            case 2    
                imagesc(reshape(x_p,Nx,Ny)) 
                %visualize result of 2D
            case 3
                plot(x_p) 
                %visuzlize the result in 3D
        end
    end
end 
%% Comparing to A\b solve
tic
A\b;
AbsB_time = toc;
%% Using Matlab's PCG function just for fun
tic
x_mat = gmres(A,b,[],tol,maxiters);
gmres_time = toc;

%% Displaying all times
disp(['TGCR: ' num2str(tgcr_time) ])
disp(['TGCR_P: ' num2str(precon_time)])
disp(['AbsB: ' num2str(AbsB_time)]) 
disp(['GMRES: ' num2str(gmres_time)])