function Jf = eval_Jf_linearSystem(x,u,p)
% evaluates the Jacobian of the vector field f() at state x, with inputs u.
% p is a structure containing all model parameters
% in particular p.A and p.B in state space model dx/dt = Ax+Bu
%
% f=eval_Jf_linear(x,u,p);

% copyright Luca Daniel, MIT 2018

Jf=p.A;
% notice this is very very simple because it is a linear function
% if your function is nonlinear you could just modify this fuction
% and make it as complicated as needed