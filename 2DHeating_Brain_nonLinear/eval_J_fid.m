function J = eval_J_fid(F,x,eps)
%This function evaluates our nonlinear Jacobian using 
%finite difference methods.  It has been written to confirm
%that we wrote our analytic Jacobian correctly.  
%~~~~~~
%Inputs
%~~~~~~
%F ~ Function handler which evaluates our forward function,
%aka, given input temperatures, compute heat flows
%x ~ N x 1, Input temperatures at which we wish to compute 
%our Jacobian
%eps ~ 1 x 1, Perturbation which we will apply to compute our 
%Jacobian
%~~~~~~~
%Outputs
%~~~~~~~
%J ~ N x N, Finite different approximation to the Jacobian

N=length(x); %number of decision variables 
J = zeros(N); %preallocate memory for or jacobian

f = F(x); %Precompute F(x) so it only needs to evaluated once

for i=1:N %loop through each input variable x_i, and compute 
    %column of the jacobian by perturbing x_i
    
    x_p = x;
    x_p(i)=x(i)+eps; %perturb x
    
    J(:,i) =(F(x_p)-f)./eps;    
end

end