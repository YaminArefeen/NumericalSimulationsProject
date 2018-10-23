function A_hat = genPrecon(flag,A)
%This function generates the appropriate inverse of our preconditioner that we want, which we can pass into TGCR_p t solve.

%Inputs
%flag ~ Indicates which preconditioner we wish to generate 
%1 => Jacobi, 2 => tridiag, 3 => lower diag
%A ~ Our M x M nodal matrix 

%Outputs
%A_hat ~ A M x M Preconditioner 

A_hat = zeros(size(A)); %preallocate space

     switch(flag)
        case 1 
               %Generate a jacobi preconditioner 
              A_hat = diag(diag(A));
             %Simply take diagonal part of nodal matrix 
         case 2
            %Generate a tridiagonal preconditioner
            N = size(A,1); %Dimension of our square matrix
            for ii = 1:N
                A_hat(ii,ii) = A(ii,ii);   
            end

            for ii = 1:N-1
                A_hat(ii,ii+1) = A(ii,ii+1);
                A_hat(ii+1,ii) = A(ii+1,ii);
            end
         case 3
             %Generate lower triangular precon
             A_hat = triu(A);
     end
end
