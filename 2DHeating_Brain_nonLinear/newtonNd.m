function [solutions,nf,ndx] = newtonNd(f,x0,c_tol,dx_tol,r_tol,t)
%this function implements a multidimensional newton algorithm
%f ~ input functio handle
%x0 ~ initial starting guess
%c_tol ~ residual convergence tolerance
%dx_tol ~ change in x convergence tolerance
%r_tol ~ relative change in x conergence tolerance
%t ~ time at which we wish to evaluate function

maxIters = 50;
nf = zeros(length(x0),maxIters); %preallocate function value at each iteration
ndx = zeros(length(x0),maxIters); %preallocate dx value at each iteration
solutions = zeros(length(x0),maxIters); %preallocate intermediate solution at each iteration

x = x0;
xprev = x;

for iter = 1:maxIters
         [func,J] = f(x,t); %evaluate our function and its jacobian 
         dx = J\-func; %solve the system for delta x
         x = x + dx;

         solutions(:,iter) = x;  
         nf(:,iter) = func; %function value at this iteration
         ndx(:,iter) = dx; %dx at this iterations

         %check for convergnece 
        if(norm(func,inf) < c_tol && norm(dx,inf) < dx_tol ... 
                && max(abs(dx))/max(abs(xprev)) < r_tol)
%                fprintf('Converged to max residual  %f in %d iterations\n',max(nf(:,iter)),iter);
               break;
        end
        xprev = x;     
end

if(iter == maxIters)
        fprintf('Nonconvergence after %d iterations\n',iter)
end

solutions = solutions(:,1:iter); 
nf = nf(:,1:iter);
ndx = ndx(:,1:iter);
%crop off all preallocated zero elements
end
