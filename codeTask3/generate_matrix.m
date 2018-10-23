function [A,b] = generate_matrix(dim_flag, visualize, N, tmsize)
%% GENERATING A NODAL MATRIX
%Dimension of the brain that we want to creat a nodal matrix fore
%1 => 1D, $2 => 2D, $3 => 3D

switch(dim_flag)
    case 1
        %Generate a 1D nodal Matrix
    case 2
        %Generate a 2D nodal Matrix
        Nx= N;%100; %Nodes in x direction
        Ny= N;%100; %Nodes in y direction
        dx = 1; %Dicretization in x direction
        dy = 1; %Dicretization in y direction
        
        % 2D Heat Process
%         tmsize = 15; %Dimension of tumor square in center
        
        [p,x_start,t_start,t_stop,max_dt_FE] = ...
            getParam_2DHeating(dx,dy,Nx,Ny,tmsize);
        
        A = p.A; %our nodal matrix
        b = p.B; %single excitation in the center
    case 3
        %Generate a 3D nodal Matrix
        N = N;%22; %Dimension in all 3, 3D dimensions
        R = zeros(3,N);
        R(1,:) = 1;
        R(2,:) = 0:N-1;
        R(3,:) = 1:N; 
        %assume all resistors are connected by value 1 in a line
        
        [E,~,alpha] = constructEIs(R,[],N); 
        G = E*alpha*E';
        %construct the 1D nodal Matrix 
        
        I = eye(N);
        G2 = sparse(kron(I,G) + kron(G,I));
        %construct the 2D nodal Matrix
        
        I2 = eye(N^2);
        A = sparse(kron(I2,G) + kron(G2,I));
        %construct the 3D Nodal Matrix 
        %for now, generate generic 3D nodal matrix just to test
        %preconditioner 
        
        b = zeros(N^3,1);
        b(floor(N^2/2),1) = 1;
        b = sparse(b);
        %apply heat to the middle of the outermost square
        %like applying heat to the edge of the brain
        
end
end