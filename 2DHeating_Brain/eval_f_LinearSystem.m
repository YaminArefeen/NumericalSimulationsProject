function f = eval_f_linearSystem(x,u,p)
% evaluates the vector field f() at state x, with inputs u.
% p is a structure containing all model parameters
% in particular p.A and p.B in state space model dx/dt = Ax+Bu
%
% f=eval_f(x,u,p);

% copyright Luca Daniel, MIT 2018

f=p.A * x + p.B * u;
