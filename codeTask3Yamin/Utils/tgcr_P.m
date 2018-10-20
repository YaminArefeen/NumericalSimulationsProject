function [x, r_norms] = tgcr_P(M,b,tol,maxiters,Ahat)
% Generalized conjugate residual method for solving Mx = b
% INPUTS
% M - matrix
% b - right hand side
% tol - convergence tolerance, terminate on norm(b - Mx) < tol * norm(b)
% maxiters - maximum number of iterations before giving up
% OUTPUTS
% x - computed solution, returns null if no convergence
% r_norms - the scaled norm of the residual at each iteration (r_norms(1) = 1)
% Ahat - Inverse of preconditioner matrix. 

[L_hat,U_hat,P_hat] = lu(sparse(Ahat));

L_hat = sparse(L_hat);
U_hat = sparse(U_hat);
P_hat = sparse(P_hat);

% Generate the initial guess for x (zero)
x = zeros(size(b));

% Set the initial residual to b - Mx^0 = b
tmp = L_hat\(P_hat*b); %Solve L (U bhat) = b for Ub^
bhat = U_hat\tmp;
r = bhat; 

% Determine the norm of the initial residual
r_norms(1) = norm(r,2);

for iter = 1:maxiters
% Use the residual as the first guess for the new
% search direction and multiply by M
  p(:,iter) = r;
  
  image = M*p(:,iter);  
  tmp = L_hat\(P_hat*image); 
  Mp(:,iter) = U_hat\tmp; 

% Make the new Mp vector orthogonal to the previous Mp vectors,
% and the p vectors M^TM orthogonal to the previous p vectors
  for j = 1:iter-1
    beta = Mp(:,iter)' * Mp(:,j);
    p(:,iter) = p(:,iter) - beta * p(:,j);
    Mp(:,iter) = Mp(:,iter) - beta * Mp(:,j);
  end;

% Make the orthogonal Mp vector of unit length, and scale the
% p vector so that M * p  is of unit length
  norm_Mp = norm(Mp(:,iter),2);
  Mp(:,iter) = Mp(:,iter)/norm_Mp;
  p(:,iter) = p(:,iter)/norm_Mp;

% Determine the optimal amount to change x in the p direction
% by projecting r onto Mp
  alpha = r' * Mp(:,iter);

% Update x and r
  x = x + alpha * p(:,iter);
  r = r - alpha * Mp(:,iter);

% Save the norm of r
  r_norms(iter+1) = norm(r,2);

% Print the norm during the iteration
% fprintf('||r||=%g i=%d\n', norms(iter+1), iter+1);

% Check convergence.
  if r_norms(iter+1) < (tol * r_norms(1))
    break;
  end
end

% Notify user of convergence
if r_norms(iter+1) > (tol * r_norms(1))
  fprintf(1, 'GCR NONCONVERGENCE!!!\n');
  x = [];
else
  fprintf(1, 'GCR converged in %d iterations\n', iter);
end

% Scale the r_norms with respect to the initial residual norm
r_norms = r_norms / r_norms(1);
